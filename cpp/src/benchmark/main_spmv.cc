
#include <cstdio>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>
#include <unordered_map>
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_spmv.hpp>

enum {KK_KERNELS, KK_KERNELS_INSP, KK_INSP};
enum {AUTO, DYNAMIC, STATIC};

typedef double Scalar;
typedef int Ordinal;
typedef int Offset;
typedef Kokkos::LayoutLeft Layout;

template<typename AType, typename XType, typename YType>
void matvec(AType& A, XType x, YType y, Ordinal rows_per_thread,
	    int team_size, int vector_length, int test, int schedule)
{
  KokkosSparse::spmv (KokkosSparse::NoTranspose, 1.0, A, x, 0.0, y);
}

int test_crs_matrix_singlevec(Ordinal numRows, Ordinal numCols,
			      int test, const char* filename,
			      Ordinal rows_per_thread, int team_size,
			      int vector_length, int schedule, int loop)
{
  typedef KokkosSparse::CrsMatrix<Scalar, Ordinal, Kokkos::DefaultExecutionSpace, void, Offset> matrix_type;
  typedef typename Kokkos::View<Scalar*, Layout> mv_type;
  typedef typename mv_type::HostMirror h_mv_type;

  srand(17312837);
  Ordinal nnzPerRow = 4;
  Scalar bandFraction = 0.001;
  std::cout << bandFraction * numRows << "\n";
  matrix_type A;
  if(filename)
    A = KokkosKernels::Impl::read_kokkos_crst_matrix<matrix_type>(filename);
  else
  {
    Offset nnz = nnzPerRow * numRows;
    A = KokkosKernels::Impl::kk_generate_sparse_matrix<matrix_type>(numRows, numCols,
								    nnz, 0, bandFraction * numCols);
  }
  numRows = A.numRows();
  numCols = A.numCols();
  Offset nnz = A.nnz();
  mv_type x("X", numCols);
  mv_type y("Y", numRows);
  h_mv_type h_x = Kokkos::create_mirror_view(x);
  h_mv_type h_y = Kokkos::create_mirror_view(y);
  h_mv_type h_y_compare = Kokkos::create_mirror(y);

  typename matrix_type::StaticCrsGraphType::HostMirror h_graph = Kokkos::create_mirror(A.graph);
  typename matrix_type::values_type::HostMirror h_values = Kokkos::create_mirror_view(A.values);

  for(int i=0; i<numCols;i++) {
    h_x(i) = (Scalar) (1.0*(rand()%40)-20.);
  }
  for(int i=0; i<numRows;i++) {
    h_y(i) = (Scalar) (1.0*(rand()%40)-20.);
  }

  // Error Check Gold Values
  for(int i=0;i<numRows;i++) {
    int start = h_graph.row_map(i);
    int end = h_graph.row_map(i+1);
    for(int j=start;j<end;j++) {
      h_values(j) = h_graph.entries(j) + i;
    }

    h_y_compare(i) = 0;
    for(int j=start;j<end;j++) {
      Scalar tmp_val = h_graph.entries(j) + i;
      int idx = h_graph.entries(j);
      h_y_compare(i)+=tmp_val*h_x(idx);
    }
  }

  Kokkos::deep_copy(x,h_x);
  Kokkos::deep_copy(y,h_y);
  Kokkos::deep_copy(A.graph.entries,h_graph.entries);
  Kokkos::deep_copy(A.values,h_values);
  mv_type x1("X1",numCols);
  Kokkos::deep_copy(x1,h_x);
  mv_type y1("Y1",numRows);

  // // //int nnz_per_row = A.nnz()/A.numRows();
  matvec(A,x1,y1,rows_per_thread,team_size,vector_length,test,schedule);

  // Error Check
  Kokkos::deep_copy(h_y,y1);
  Scalar error = 0;
  Scalar sum = 0;
  for(int i=0;i<numRows;i++) {

    error += (h_y_compare(i)-h_y(i))*(h_y_compare(i)-h_y(i));
    sum += h_y_compare(i)*h_y_compare(i);
  }

  int num_errors = 0;
  double total_error = 0;
  double total_sum = 0;
  num_errors += (error/(sum==0?1:sum))>1e-5?1:0;
  total_error += error;
  total_sum += sum;

  // Benchmark
  double min_time = 1.0e32;
  double max_time = 0.0;
  double ave_time = 0.0;
  for(int i=0;i<loop;i++) {
    Kokkos::Timer timer;
    matvec(A,x1,y1,rows_per_thread,team_size,vector_length,test,schedule);
    Kokkos::fence();
    double time = timer.seconds();
    ave_time += time;
    if(time>max_time) max_time = time;
    if(time<min_time) min_time = time;
  }

  // Performance Output
  double matrix_size = 1.0*((nnz*(sizeof(Scalar)+sizeof(Ordinal)) + numRows*sizeof(Offset)))/1024/1024;
  double vector_size = 2.0*numRows*sizeof(Scalar)/1024/1024;
  double vector_readwrite = (nnz+numRows)*sizeof(Scalar)/1024/1024;

  double problem_size = matrix_size+vector_size;
  printf("NNZ NumRows NumCols ProblemSize(MB) AveBandwidth(GB/s) MinBandwidth(GB/s) MaxBandwidth(GB/s) AveGFlop MinGFlop MaxGFlop aveTime(ms) maxTime(ms) minTime(ms) numErrors\n");
  printf("%i %i %i %6.2lf ( %6.2lf %6.2lf %6.2lf ) ( %6.3lf %6.3lf %6.3lf ) ( %6.3lf %6.3lf %6.3lf ) %i RESULT\n",nnz, numRows,numCols,problem_size,
          (matrix_size+vector_readwrite)/ave_time*loop/1024, (matrix_size+vector_readwrite)/max_time/1024,(matrix_size+vector_readwrite)/min_time/1024,
          2.0*nnz*loop/ave_time/1e9, 2.0*nnz/max_time/1e9, 2.0*nnz/min_time/1e9,
          ave_time/loop*1000, max_time*1000, min_time*1000,
          num_errors);
  return (int)total_error;
}

void print_help() {
  printf("SPMV benchmark code written by Christian Trott.\n");
  printf("OpenMP implementations written by Simon Hammond (Sandia National Laboratories).\n\n");
  printf("Options:\n");
  printf("  --rows [N]  : generate a semi-random banded (band size 0.01xN) NxN matrix\n");
  printf("  --cols [N]  : generate a semi-random banded (band size 0.01xN) NxN matrix\n");
  printf("  --test [OPTION] : Use different kernel implementations\n");
  printf("                    Options:\n");
  printf("                      kk,kk-kernels          (Kokkos/Trilinos)\n");
  printf("                      kk-insp                (Kokkos Structure Inspection)\n");
  printf("  --schedule [SCH]: Set schedule for kk variant (static,dynamic,auto [ default ]).\n");
  printf("  -f [file]       : Read in Matrix Market formatted text file 'file'.\n");
  printf("  -fb [file]      : Read in binary Matrix files 'file'.\n");
  printf("  --write-binary  : In combination with -f, generate binary files.\n");
  printf("  -rpt [K]        : Number of Rows assigned to a thread.\n");
  printf("  -ts [T]         : Number of threads per team.\n");
  printf("  -vl [V]         : Vector-length (i.e. how many Cuda threads are a Kokkos 'thread').\n");
  printf("  -l [LOOP]       : How many spmv to run to aggregate average time. \n");
}

int main(int argc, char **argv)
{
 long long int rows = 110503; // a prime number
 long long int cols = 110503; // a prime number
 //int numVecs = 4;
 int test=KK_KERNELS;
 //int type=-1;
 char* filename = NULL;

 int rows_per_thread = -1;
 int vector_length = -1;
 int team_size = -1;
 int schedule=AUTO;
 int loop = 100;

 if(argc == 1) {
   print_help();
   return 0;
 }

 for(int i=0;i<argc;i++)
 {
   //if((strcmp(argv[i],"-s")==0)) {size=atoi(argv[++i]); continue;}
   if((strcmp(argv[i],"--rows")==0)) {rows=atoi(argv[++i]); continue;}
   if((strcmp(argv[i],"--cols")==0)) {cols=atoi(argv[++i]); continue;}

   //if((strcmp(argv[i],"-v")==0)) {numVecs=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--test")==0)) {
    i++;
    if(i == argc)
    {
      std::cerr << "Must pass algorithm name after '--test'";
      exit(1);
    }
    if((strcmp(argv[i],"kk-kernels")==0))
      test = KK_KERNELS;
    if((strcmp(argv[i],"kk-kernels-insp")==0))
      test = KK_KERNELS_INSP;
    if((strcmp(argv[i],"kk-insp")==0))
      test = KK_INSP;
    continue;
  }
  //if((strcmp(argv[i],"--type")==0)) {type=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-f")==0)) {filename = argv[++i]; continue;}
  if((strcmp(argv[i],"-fb")==0)) {filename = argv[++i]; continue;}
  if((strcmp(argv[i],"-rpt")==0)) {rows_per_thread=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-ts")==0)) {team_size=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-vl")==0)) {vector_length=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-l")==0)) {loop=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--schedule")==0)) {
    i++;
    if((strcmp(argv[i],"auto")==0))
      schedule = AUTO;
    if((strcmp(argv[i],"dynamic")==0))
      schedule = DYNAMIC;
    if((strcmp(argv[i],"static")==0))
      schedule = STATIC;
    continue;
  }
  if((strcmp(argv[i],"--help")==0) || (strcmp(argv[i],"-h")==0)) {
    print_help();
    return 0;
  }
 }

 Kokkos::initialize(argc,argv);

 int total_errors = test_crs_matrix_singlevec(rows,cols,test,filename,
					      rows_per_thread,
					      team_size,vector_length,schedule,loop);

 if(total_errors == 0)
   printf("Kokkos::MultiVector Test: Passed\n");
 else
   printf("Kokkos::MultiVector Test: Failed\n");


  Kokkos::finalize();
}
