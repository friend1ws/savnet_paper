letterA <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  x <- c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6)
  y <- c(0,10,10,0,0,3,3,0,0,4,7.5,4,4)
  x <- 0.1*x
  y <- 0.1*y
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  if (is.null(id)){
    id <- c(rep(1,9),rep(2,4))
    fill <- c(rep("A", 9),rep("white", 4))
  }else{
    id <- c(rep(id,9),rep(id+1,4))
    fill <- c(rep("A", 9),rep("white", 4))
  }
  
  data.frame(x=x,y=y,id=id,fill=fill)
}

## T
letterT <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  x <- c(0,10,10,6,6,4,4,0)
  y <- c(10,10,9,9,0,0,9,9)
  x <- 0.1*x
  y <- 0.1*y
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  if (is.null(id)){
    id <- rep(1,8)
    fill <- rep("T", 8)
  }else{
    id <- rep(id,8)
    fill <- rep("T", 8)
  }
  
  data.frame(x=x,y=y,id=id,fill=fill)
}

## C
letterC <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  if (is.null(id)){
    id <- rep(1,length(x))
    fill <- rep("C", length(x))
  }else{
    id <- rep(id,length(x))
    fill <- rep("C", length(x))
  }
  
  data.frame(x=x,y=y,id=id,fill=fill)
}


## G
letterG <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  
  h1 <- max(y.l1)
  r1 <- max(x.l1)
  
  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)
  
  
  
  if (is.null(id)){
    id <- c(rep(1,length(x)),rep(2,length(x.add)))
    fill <- rep("G", length(x) + length(x.add))
  }else{
    id <- c(rep(id,length(x)),rep(id+1,length(x.add)))
    fill <- rep("G", length(x) + length(x.add))
  }
  
  x <- c(rev(x),x.add)
  y <- c(rev(y),y.add)
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  data.frame(x=x,y=y,id=id,fill=fill)
  
}


M <- matrix(rgamma(4 * 8, 0.1, 1), 4, 8)


gg_seq_logo <- function(M, is.scale = FALSE, A_col = "#00BA38", C_col = "#619CFF", G_col = "#B79F00", T_col = "#F8766D") {
 
  M <- sweep(M, 2, colSums(M), FUN="/")
  
  if (is.scale) {
    for (i in 1:ncol(M)) {
      ic <- 2 + sum(sapply(M[, i], function(x) {ifelse(x > 0, x*log2(x), 0)}))
      M[, i] <- M[, i] * ic
    } 
  }
  
  letter_data <- c()
  for(i in 1:ncol(M)) {
    cur_y <- 0
    cur_id <- 8 * (i - 1) 
  
    letter_order <- order(M[,i])
    for(tl in letter_order) {
      if (tl == 1) {
        letter_data <- rbind(letter_data, letterA(i - 0.49, cur_y, 0.97 * M[1, i], 0.98, cur_id))
        cur_y <- cur_y + 0.97 * M[1, i] + 0.01
      } else if (tl == 2) {
        letter_data <- rbind(letter_data, letterC(i - 0.49, cur_y, 0.97 * M[2, i], 0.98, cur_id + 2))
        cur_y <- cur_y + 0.97 * M[2, i] + 0.01   
      } else if (tl == 3) {
        letter_data <- rbind(letter_data, letterG(i - 0.49, cur_y, 0.97 * M[3, i], 0.98, cur_id + 4))
        cur_y <- cur_y + 0.97 * M[3, i] + 0.01         
      } else if (tl == 4) {
        letter_data <- rbind(letter_data, letterT(i - 0.49, cur_y, 0.97 * M[4, i], 0.98, cur_id  +6))
        cur_y <- cur_y + 0.97 * M[4, i] + 0.01         
      } else {
        stop("tl must be any of 1 to 4")
      }
    }
  }

  ggplot(letter_data, aes(x = x, y = y, group = id, fill = fill)) + geom_polygon() +
    scale_fill_manual(values = c("A" = A_col, "C" = C_col, "G" = G_col, "T" = T_col, "white" = "#FFFFFF")) +
    guides(fill = FALSE) +
    labs(x = "Position", y = "Probability")

}



