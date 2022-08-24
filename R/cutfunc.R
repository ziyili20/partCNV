cutfunc <- function(choice, lower = 0, upper = 1) {
    newchoice <- choice
    newchoice[choice<=lower] <- 0.01
    newchoice[choice>=upper] <- 0.99
    return(newchoice)
}
