#!/bin/bash
import json

# load the version number from package.json
version = json.load(open('package.json'))['version']
print('* found version number is %s' % version)

# filenames
fn_src = './src/meta.js'
fn_out_esm = './dist/metajs-esm-%s.js' % version
fn_out_esm_latest = './dist/metajs-esm-latest.js'
fn_out_script = './dist/metajs-%s.js' % version
fn_out_script_latest = './dist/metajs-latest.js'

# read the src
raw = open(fn_src).read()

###########################################################
# convert the ESM version
###########################################################

# replace the content
# import { create, all } from 'mathjs'
# const math = create(all, {});
txt = raw.replace("import { create, all } from 'mathjs'", "")
txt = txt.replace("const math = create(all, {});", "")

# output esm
open(fn_out_esm, 'w').write(txt)
open(fn_out_esm_latest, 'w').write(txt)

print('* output %s -> %s + %s' % (
    fn_src,
    fn_out_esm,
    fn_out_esm_latest
))


###########################################################
# convert the script version
###########################################################

txt = txt.replace("'use strict';", "// Meta.js")
txt = txt.replace("export const ", "var ")

# output esm
open(fn_out_script, 'w').write(txt)
open(fn_out_script_latest, 'w').write(txt)

print('* output %s -> %s + %s' % (
    fn_src,
    fn_out_script,
    fn_out_script_latest
))