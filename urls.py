from django.conf.urls import patterns, include, url
import settings
from django.contrib import admin
from primer3.primerviews import *
from django.contrib.staticfiles.urls import staticfiles_urlpatterns  
admin.autodiscover()

urlpatterns = patterns('',
    url(r'^site_medias/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT }),   
    url(r'^static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.STATIC_ROOT }),
    url(r'^admin/', include(admin.site.urls)),
    url(r'^primer/',primerdesigner),
    url(r'^seqview/',seqview),
    url(r'^basemodify/',basemodify),
    url(r'^seqinput/',seqinput),
    url(r'^seqread/',seqread),
    url(r'^pcr-principle/',pcrprincp),
    url(r'^primer-manual/',manual),
    url(r'^multipcr/',multiprimer),
    url(r'^multiseqview/',multiseqview),
    
)
