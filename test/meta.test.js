'use strict';

import { metajs } from '../src/meta.js';
import assert from 'assert'

describe('testing metajs base functions', () => {
    it('backtransf 1 OR null null should be 1', () => {

        assert.equal(
            metajs.backtransf(1, 'OR', null, null), 
            1
        );
    });

    // related to ci95
    it('ci95(0.5, 0.2) should be ', () => {
        assert.deepEqual(
            metajs.ci95(0.5, 0.2), 
            {
                "TE": 0.5,
                "seTE": 0.2,
                "pval": 0.01242,
                "level": 0.95,
                "lower": 0.108,
                "upper": 0.892,
                "statistic": 2.5
            } 
        );
    });

    // just an empty test
    it('at least one pass', () => {
        assert.equal(
            1,
            1
        );
    });
});
