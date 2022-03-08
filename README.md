# Meta.js

Meta.js is designed to be a JavaScript library for online meta-analysis.
It's built on [math.js](https://mathjs.org/) and [NumJs](https://github.com/nicolaspanel/numjs) to support calculation tasks related meta-analysis.
Meta.js is still in very early development stage, please don't use it for production purpose.

Since Meta.js is purely built on JavaScript, it doesn't require any runtime installation, which means that you don't need to install R, Python, or any other libraries.
Meta.js can be used to build serverless web applications, and we present a prototype to demonstrate how it looks:

[https://ohnlp.github.io/Meta.js/](https://ohnlp.github.io/Meta.js/)

![Meta.js Screenshot](https://raw.githubusercontent.com/OHNLP/Meta.js/main/static/img/screenshot.png)

## Methods

At present, Meta.js implements the basic functions for meta-analysis of single proportions (metaprop), and meta-analysis of binary outcome data (metabin), but it only supports very limited methods and measures.

## Performance

To validate the correctness and evaluate the performance of Meta.js, the following files are included:

1. `make_sample.py`: this Python script can make test dataset of binary outcomes for test purpose.
2. `pwma.R`: this R script provides the performance test. We test it in R 3.6.3 (x86_64-pc-linux-gnu) with `meta` 4.18-0 version.
3. `sample.csv`: the sample data file for test, which contains 10,000 outcomes, and each outcome contains random number of studies (1-100)

