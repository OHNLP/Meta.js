'use strict';

import { metajs } from './src/meta.js';
import Papa from 'papaparse';
import fs from 'fs';
import * as dfd from "danfojs-node";

console.log(
    metajs.expit(1)
)

// locate the input stream
let f = fs.createReadStream('./testsets/test_input.csv');
let fo = fs.createReadStream('./testsets/test_result_OR.csv');

var rs1 = null;
var rs2 = null;

Papa.parse(f, {
    header: true,
    comments: "#",
    complete: function(rs) {
        console.log('* loaded ' + rs.data.length + ' records in input testset');
        rs1 = rs;
        Papa.parse(fo, {
            header: true,
            comments: "#",
            complete: function(rs) {
                console.log('* loaded ' + rs.data.length + ' records in result testset');
                rs2 = rs;

                // ok, now let's test
                do_something()
            }
        });
    }
});

function do_something() {
    // create a dataframe
   var df = new dfd.DataFrame(rs1.data);
   var oc_names = df['outcome'].unique().values;

   // get the dataframe for results
   var dfr = new dfd.DataFrame(rs2.data);

   dfr.print();
}