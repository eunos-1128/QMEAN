"""website URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.8/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Add an import:  from blog import urls as blog_urls
    2. Add a URL to urlpatterns:  url(r'^blog/', include(blog_urls))
"""
from django.conf.urls import include, url
from django.contrib import admin
from django.views.generic.base import TemplateView
import views

urlpatterns = [
	url(r'^/?$', views.home, name='home'),
	url(r'^structure_upload/$', views.structure_upload, name='structure_upload'),
	url(r'^sequence_upload/$', views.sequence_upload, name='sequence_upload'),
    url(r'^uploaded_structure/(?P<projectid>\w{6})/(?P<modelid>model_[\d]{3}).pdb$', views.uploaded_structure, name='uploaded_structure'),
    url(r'^model_structure/(?P<projectid>\w{6})/(?P<modelid>model_[\d]{3}).pdb$', views.model_structure, name='model_structure'),
    url(r'^results/(?P<projectid>\w{6}).zip$', views.archive, name='archive'),
    url(r'^results/(?P<projectid>\w{6})/(?P<modelid>model_[\d]{3}).zip$', views.archive, name='archive'),
	url(r'^results/(?P<projectid>\w{6})/$', views.results, name='results'),
    url(r'^project_status/(?P<projectid>\w{6})/$', views.get_project_status_json, name='project_status'),
    url(r'^qmean_pix/(?P<projectid>\w{6})/(?P<modelid>model_[\d]{3})/(?P<name>[\w\._]+)$', views.qmean_pix, name='qmean_pix'),
    url(r'^help/?', views.help, name='help'),
    url(r'^references/?', views.references, name='references'),
    url(r'^contact/?', views.contact, name='contact'),
    url(r'^404/$',  TemplateView.as_view(template_name='404.html')),
    url(r'^500/$',  TemplateView.as_view(template_name='500.html')),
]
