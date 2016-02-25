## Run the tests if the system variable 'R_CHECK_TESTS' is set to TRUE

# code from bioconductor package 'derfinder' -- thank you!

flag <- as.logical(Sys.getenv('R_CHECK_TESTS'))
if(!is.na(flag)) {
    if(flag == TRUE) {
        library('testthat')
        test_check('test_ebm_main')
    }
} else {
    message('Set the system variable R_CHECK_TESTS to TRUE to run all the tests.')
}
