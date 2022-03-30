#!/usr/bin/env python

# Copyright (C) 2019 Thermo Fisher Scientific. All Rights Reserved

# vcMerge plugin on Genexus
# 2022/3/29

import glob
import sys
import subprocess
import json
import os
import re
import shutil
from ion.plugin import *
from django.utils.functional import cached_property
from django.conf import settings
from django.template.loader import render_to_string




class vcMerge_GX(IonPlugin):
  version = '1.0.0.0'
  author = "longfei.fu@thermofisher.com"
  runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]

  # a simple cached version of the start plugin property
  @cached_property
  def startplugin_json(self):
    return self.startplugin

  def merge(self):
    plugin_result_dir = self.startplugin_json['runinfo']['plugin'].get('results_dir')
    lane_info_dir     = self.startplugin_json['runinfo'].get('analysis_dir')
    results_name      = self.startplugin_json['expmeta'].get('results_name')
    outfile = "%s/%s.TSVC_variants.merged.vcf.xls" % (plugin_result_dir,results_name)

    print "plugin result dir is: %s" % (plugin_result_dir)
    print "lane dir is: %s" % (lane_info_dir)
    print "outfile is: %s" % (outfile)

    abs_path = os.path.abspath(__file__)
    this_dir = os.path.dirname(abs_path)
    #cmd = "perl %s/merge_main.pl %s %s" % (this_dir,lane_info_dir,plugin_result_dir)
    cmd = "perl %s/merge_main.pl %s %s" % (this_dir,lane_info_dir,outfile)
    print "cmd is %s" % (cmd)

    print "Start running the vcMerge_GX plugin."
    os.system(cmd)
    print "Finished the vcMerge_GX plugin."

  def launch(self,data=None):
    self.merge()
    with open("vcMerge_GX_block.html","w") as f:
      f.write('<html><body>To download: "Right Click" -> "Save Link As..."<br>\n')
      for merge_vcf in glob.glob('*.merged.vcf.xls'):
        print(merge_vcf)
        f.write('<a href="%s">%s</a><br>\n'
                % (merge_vcf,merge_vcf))
        f.write('</body></html>\n')
    
    return True


if __name__ == "__main__":
    PluginCLI()
