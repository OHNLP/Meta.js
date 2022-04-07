#!/bin/bash

version = '0.0.3'
fn_src = './src/meta.js'
fn_out_version = './dist/metajs-esm-%s.js' % version
fn_out_latest = './dist/metajs-esm-latest.js'

# read the src
txt = open(fn_src).read()

# replace the content
# import { create, all } from 'mathjs'
# const math = create(all, {});
txt = txt.replace("import { create, all } from 'mathjs'", "")
txt = txt.replace("const math = create(all, {});", "")

# output
open(fn_out_version, 'w').write(txt)
open(fn_out_latest, 'w').write(txt)

print('* output %s -> %s + %s' % (
    fn_src,
    fn_out_version,
    fn_out_latest
))