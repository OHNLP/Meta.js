'use strict';

import { pchisq } from '../src/meta.js';
import assert from 'assert'


describe('testing pchisq function', () => {
    it('pchisq.compute(0, 2, false) should be 1', () => {
        assert.equal(
            pchisq.compute(0, 2, false),
            1
        );
    });
    it('pchisq.compute(0.3, 2, false) should be 0.86071', () => {
        assert.equal(
            pchisq.compute(0.3, 2, false),
            0.86071
        );
    });
    it('pchisq.compute(0.3, 3, false) should be 0.96003', () => {
        assert.equal(
            pchisq.compute(0.3, 3, false),
            0.96003
        );
    });
    it('pchisq.compute(0.9, 5, false) should be 0.97022', () => {
        assert.equal(
            pchisq.compute(0.9, 5, false),
            0.97022
        );
    });
});