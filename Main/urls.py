from django.urls import path, include

from django.contrib import admin

admin.autodiscover()

import Core.views

# To add a new path, first import the app:
# import blog
#
# Then add the new path:
# path('blog/', blog.urls, name="blog")
#
# Learn more here: https://docs.djangoproject.com/en/2.1/topics/http/urls/


urlpatterns = [
    path("admin/", admin.site.urls),
    path("", Core.views.home),
    path("snp/", Core.views.snp),
    path("primer/", Core.views.primerdesighn),
    path("chemical/", Core.views.chemical),
    path("genes/", Core.views.genes),
    path("table/", Core.views.table),
    path("phy/", Core.views.phylo),
    path("phygene/", Core.views.phylogene),
    path("pathway/", Core.views.pathway),
    path("disease/", Core.views.disease_database),
    path("dcatcher/", Core.views.disease_catcher),
    path("srange/", Core.views.snp_range),

]
