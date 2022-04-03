'use strict';

import { pnorm } from '../src/meta.js';
import assert from 'assert'

describe('testing pnorm function', () => {
    it('pnorm.compute(0, false) should be 0.5', () => {
        assert.equal(
            pnorm.compute(0, false),
            0.5
        );
    });
    it('pnorm.compute(1, false) should be 0.15866', () => {
        assert.equal(
            pnorm.compute(1, false),
            0.15866
        );
    });
    it('pnorm.compute(2, false) should be 0.02275', () => {
        assert.equal(
            pnorm.compute(2, false),
            0.02275
        );
    });
});
