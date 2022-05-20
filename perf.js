'use strict'
import { performance } from 'node:perf_hooks';
import { metajs } from './src/meta.js';

function simple_test() {
    var rs = [
        {study: 'SA', sm: 0.62, lower: 0.46, upper: 0.82, t1: 'CaboNivo', t2: 'Suni', year: 2020},
        {study: 'SB', sm: 0.54, lower: 0.46, upper: 0.63, t1: 'NivoIpi', t2: 'Suni', year: 2021},
        {study: 'SC', sm: 0.68, lower: 0.58, upper: 0.81, t1: 'AteBev', t2: 'Suni', year: 2022},
        {study: 'SD', sm: 0.80, lower: 0.62, upper: 1.03, t1: 'Pazo', t2: 'Suni', year: 2022},
        {study: 'SE', sm: 0.74, lower: 0.45, upper: 1.20, t1: 'AteBev', t2: 'Suni', year: 2022},
        {study: 'SF', sm: 1.12, lower: 0.91, upper: 1.38, t1: 'PemAxi', t2: 'Suni', year: 2022},
    ];

    performance.mark("example-start");
    console.log('* started!')
    const t0 = performance.now();

    var n_times = 1000;
    for (let i = 0; i < n_times; i++) {
        var ret = metajs.netmeta(rs, {});
    }
    const t1 = performance.now();
    const dur = (t1 - t0) / 1000;
    console.log('* ' + n_times + ': ' + dur);
    
    performance.mark("example-end");
    performance.measure("example", "example-start", "example-end");
    console.log('* done');
}

simple_test();