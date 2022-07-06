# Meta.js

Meta.js is designed to be a JavaScript library for online real-time meta-analysis.
It's built on [math.js](https://mathjs.org/) to support calculation tasks related meta-analysis.
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
2. `pwma.R`: this R script provides the performance test.
3. `sample.csv`: the sample data file for test, which contains 10,000 outcomes, and each outcome contains random number of studies (1-100)
4. `gen_test_output.r`: this R script provides the ground truth for testing.

To evaluate the performance of Meta.js, we conducted two experiments, 1) calculation of the fixed effects on randomly selected 1,000 outcomes for 50 times, and 2) calculation of the fixed effects of different number of outcomes. In both R and JavaScript versions, the test dataset is pre-loaded, and the calculation time is recorded.

We tested our program on a machine with Intel Core i5-10400 and 16G ram, running Ubuntu . The Google Chrome is Version 99.0.4844.51 (Official 64-bit). The R environment is 3.6.3 with `meta` library 4.18-0 and RStudio 1.4.1106. 

The preliminary results are as follows:

<table style="width:100%;">
<tr>
  <td><img src="https://raw.githubusercontent.com/OHNLP/Meta.js/main/static/img/perf-RvsJS.png"></td>
  <td><img src="https://raw.githubusercontent.com/OHNLP/Meta.js/main/static/img/perf-10k.png"></td>
</tr>
<tr>
  <td colspan="2">* the actual value may vary on differnet machines</td>
</tr>
</table>

The median time of the JavaScript version is 1.12s while the R version is 17.94s for 1000 outcomes. As the number of outcomes increases, the absolute differences in time also increase.

# License

Apache-2.0 License

# Usage

Before loading Meta.js, please ensure the math.js is loaded.
Meta.js can be loaded as a regular JavaScript file in the browser, use the global variable math to access the libary once loaded:

```
<!DOCTYPE HTML>
<html>
<head>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/mathjs/10.1.1/math.min.js" type="text/javascript"></script>
  <script src="meta.js" type="text/javascript"></script>
</head>
<body>
  <script type="text/javascript">
    var rs = [
        [12,393,2,396, 'S1', 'T','C'],
        [24,230,24,281, 'S2', 'T','C'],
    ]
    
    var rst = metajs.metabin(
        rs,
        {
            'sm': 'OR'
        }
    );
    
    console.log(rst);
  </script>
</body>
</html>
```

# Development

`Meta.js` leverages Node.js and [Mocha](https://mochajs.org/) to conduct unittests to ensure the correctness of calculation. 
So you need to ensure a latest Node.js is installed before development. 
Our development environment is Ubuntu 20.04.
You can also install Ubuntu 20.04 VM on Windows to ensure the program works.

After clone your own copy, run the following command to install the dependency:

```
npm install
```

Then, you can run `npm test` to check the test results.

# Change Log

## 0.0.4 (2022-04-14)

- Added checking for both zero cases
- Added metaprop for incidence analysis
- Added make script for creating ESM and script tag export
- Updated test cases
- Updated demo page with more columns
- Fixed mis-alignments for expanding results

## 0.0.3 (2022-03-31)

- Added metabin with random effect
- Added random effect estimation by DL method
- Updated dependency
- Updated using Node.js for unit tests

## 0.0.2

- Added support for basic math packages

## 0.0.1

- Added basic functions for fixed effect estimation