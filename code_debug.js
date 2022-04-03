'use strict';

import { metajs } from './src/meta.js';

var rs = [
    [0,393,2,396, 'S1', 'T','C'],
    [24,230,24,281, 'S2', 'T','C'],
]

var rst = metajs.metabin(rs);

console.log(rst);