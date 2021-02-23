import re
import shutil
import logging

AUTHOR = 'Francesco  Rizzi'

# M_SITE_LOGO = 'gr/logo.png'
M_SITE_LOGO_TEXT = 'SHWAV'

SITENAME = 'ShWav'
SITESUBTITLE = 'Elastic Shear Waves'
# SITEURL = 'http://francescorizzi.net/SHAW/output'
SITEURL = 'https://github.com/fnrizzi/SHAW/docs'

# M_BLOG_NAME = ''
# M_BLOG_URL = 'blog/'

PATH = 'content'

STATIC_URL = 'static/{path}'
STATIC_SAVE_AS = 'static/{path}'
STATIC_PATHS = ['img', 'showcase']
EXTRA_PATH_METADATA = {'img/favicon.ico': {'path': '../favicon.ico'}}

# ARTICLE_PATHS = ['blog']
# ARTICLE_EXCLUDES = ['blog/authors', 'blog/categories', 'blog/tags']

PAGE_PATHS = ['']
PAGE_EXCLUDES = ['doc', 'img']
READERS = {'html': None} # HTML files are only ever included from reST

PAGE_URL = '{slug}/'
PAGE_SAVE_AS = '{slug}/index.html'

# ARCHIVES_URL = 'blog/'
# ARCHIVES_SAVE_AS = 'blog/index.html'
# ARTICLE_URL = 'blog/{slug}/' # category/ is part of the slug
# ARTICLE_SAVE_AS = 'blog/{slug}/index.html'
# DRAFT_URL = 'blog/{slug}/' # so the URL is the final one
# DRAFT_SAVE_AS = 'blog/{slug}/index.html'
# AUTHOR_URL = 'blog/author/{slug}/'
# AUTHOR_SAVE_AS = 'blog/author/{slug}/index.html'
# CATEGORY_URL = '{slug}/'
# CATEGORY_SAVE_AS = '' #'{slug}/index.html'
# TAG_URL = 'tag/{slug}/'
# TAG_SAVE_AS = '' #'blog/tag/{slug}/index.html'

AUTHORS_SAVE_AS = None # Not used
CATEGORIES_SAVE_AS = None # Not used
TAGS_SAVE_AS = None # Not used

PAGINATION_PATTERNS = [(1, '{base_name}/', '{base_name}/index.html'),
                       (2, '{base_name}/{number}/', '{base_name}/{number}/index.html')]

TIMEZONE = 'Europe/Rome'

DEFAULT_LANG = 'en'

import platform
if platform.system() == 'Windows':
    DATE_FORMATS = {'en': ('usa', '%b %d, %Y')}
else:
    DATE_FORMATS = {'en': ('en_US.UTF-8', '%b %d, %Y')}

FEED_ATOM = None
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

M_LINKS_NAVBAR1 = [('Build', '', '',
                    [
                      ("Host Serial Kokkos Build", 'build/kokkos_host_serial', 'build/kokkos_host_serial'),
                      ("Host OpenMP Kokkos Build", 'build/kokkos_host_omp',    'build/kokkos_host_omp')
                    ]),

                    ('Get Started', '', '',
                    [
                      ("Equations and discretization", 'getstarted/goveq',          'getstarted/goveq'),
                      ('Input File',                   'getstarted/inputfile',      'getstarted/inputfile'),
                      ('Material Models',              'getstarted/materialmodels', 'getstarted/materialmodels'),
                    ]),

                   ('Demos', '', '',
                    [
                      ("Single Forcing Run",            'demos/rank1fom',      'demos/rank1fom'),
                      ("Multi-forcing Run with rank-1", 'demos/rank1fomMulti', 'demos/rank1fomMulti'),
                      ("Multi-forcing Run wirh rank-2", 'demos/rank2fom',      '')
                    ]),

                   ('Various', '', '', [("License", 'license', 'license')])]
M_LINKS_NAVBAR2 = []
M_LINKS_FOOTER1 = []
M_LINKS_FOOTER2 = []
M_LINKS_FOOTER3 = []
M_LINKS_FOOTER4 = []


M_CSS_FILES = ['https://fonts.googleapis.com/css?family=Source+Code+Pro:400,400i,600%7CSource+Sans+Pro:400,400i,600,600i&subset=latin-ext',
               'static/m-dark.css']

M_FINE_PRINT = """
| Site created by Francesco Rizzi. Powered by `Pelican <https://getpelican.com>`_
  and `m.css <https://mcss.mosra.cz>`_, using a theme adapted from `magnum <https://magnum.graphics/>`_.
| This site content is `available on GitHub <https://github.com/fnrizzi/ElasticShearWaves>`_.
"""

# M_NEWS_ON_INDEX = ("Latest news on our blog", 3)

DEFAULT_PAGINATION = 10

PLUGIN_PATHS = ['m.css/plugins']
PLUGINS = ['m.abbr',
           'm.alias',
           'm.code',
           'm.components',
           'm.dot',
           'm.dox',
           'm.filesize',
           'm.gh',
           'm.gl',
           'm.htmlsanity',
           'm.images',
           'm.link',
           'm.math',
           'm.metadata',
           'm.plots',
           'm.sphinx',
           'm.vk']

FORMATTED_FIELDS = ['summary', 'description', 'landing', 'badge', 'header', 'footer']

THEME = 'm.css/pelican-theme/'
THEME_STATIC_DIR = 'static/'

M_THEME_COLOR = '#22272e'
# M_SOCIAL_TWITTER_SITE = ''
# M_SOCIAL_TWITTER_SITE_ID =
# M_SOCIAL_IMAGE = ''
# M_SHOW_AUTHOR_LIST = True

M_HTMLSANITY_SMART_QUOTES = True
M_HTMLSANITY_HYPHENATION = True

_magnum_colors_src = re.compile(r"""<span class="mh">0x(?P<hex>[0-9a-f]{6})(?P<alpha>[0-9a-f]{2})?(?P<literal>_s?rgba?f?)</span>""")
_magnum_colors_dst = r"""<span class="mh">0x\g<hex>\g<alpha>\g<literal><span class="m-code-color" style="background-color: #\g<hex>;"></span></span>"""
# #_zwnj_in_console_colors_src = re.compile(

M_CODE_FILTERS_POST = {
    'C++': lambda str: _magnum_colors_src.sub(_magnum_colors_dst, str)
}

if not shutil.which('latex'):
    logging.warning("LaTeX not found, fallback to rendering math as code")
    M_MATH_RENDER_AS_CODE = True

DIRECT_TEMPLATES = ['archives']
PAGINATED_TEMPLATES = {'archives': None, 'tag': None, 'category': None, 'author': None}

SLUGIFY_SOURCE = 'basename'
# PATH_METADATA = '(?P<slug>.+).rst'
PATH_METADATA = '(blog/)?(?P<slug>.+).rst'
# PATH_METADATA = '(?P<slug>.*)\..*'
SLUG_REGEX_SUBSTITUTIONS = [
        (r'[^\w\s-]', ''),  # remove non-alphabetical/whitespace/'-' chars
        (r'(?u)\A\s*', ''),  # strip leading whitespace
        (r'(?u)\s*\Z', ''),  # strip trailing whitespace
        (r'[-\s]+', '-'),  # reduce multiple whitespace or '-' to single '-'
        (r'C\+\+', 'cpp'),
    ]
