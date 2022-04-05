module.exports = {
    entry: "./src/meta.js",
    output: {
        filename: "./metajs-0.0.2.js",
        library: {
            name: 'metajs',
            type: 'umd'
        }
    },
    mode: 'development',
    // Enable sourcemaps for debugging webpack's output.
    devtool: "source-map",
};