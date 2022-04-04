'use strict';

import { metajs } from './src/meta.js';

function test_metabin() {
    var rs = [
        [0,393,0,396, 'S1', 'T','C'],
        // [24,230,24,281, 'S2', 'T','C'],
    ]
    
    var rst = metajs.metabin(
        rs,
        {
            'sm': 'OR'
        }
    );
    
    console.log(rst);
}

function test_metaprop() {
    var rs = [
        [2, 20,  'S1'],
        [5, 90,  'S2'],
        [20,100, 'S3'],
    ];

    var rst = metajs.metaprop(
        rs, 
        {
            'sm': 'PFT'
        }
    );

    console.log(rst);
}

test_metaprop();