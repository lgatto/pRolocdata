test_that("data validity", {
    items <- pRolocdata()$results[, "Item"]
    e <- environment()
    for (it in items) {
        data(list = it, envir = e)
        expect_true(validObject(get(it, envir = e)))
    }
})
