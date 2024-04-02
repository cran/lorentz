### R code from vignette source 'lorentz.Rnw'

###################################################
### code chunk number 1: lorentz.Rnw:90-92
###################################################
library("lorentz")
library("magrittr")


###################################################
### code chunk number 2: floater
###################################################
                          # NB: speed of light = 1 by default
u <- as.3vel(c(0.6,0,0))  # coerce to a three-velocity
u
as.4vel(u)                # four-velocity is better for calculations
(B <- boost(u))           # transformation matrix


###################################################
### code chunk number 3: flighter
###################################################
B %*% c(1,0,0,0)


###################################################
### code chunk number 4: flayter
###################################################
B %*% c(4,6,-8,9)


###################################################
### code chunk number 5: flooter
###################################################
B %*% c(1,1,0,0)
B %*% c(1,-1,0,0)


###################################################
### code chunk number 6: lorentz.Rnw:297-299
###################################################
(B <- boost(as.3vel(c(0.8,0,0))))  # 0.8c left to right
solve(B) %*% c(1,0,0,0)            # active transform


###################################################
### code chunk number 7: lorentz.Rnw:308-312
###################################################
u <- as.3vel(c(0.3,-0.4,+0.8))
v <- as.3vel(c(0.4,+0.2,-0.1))
L <- boost(u) %*% boost(v)
L


###################################################
### code chunk number 8: lorentz.Rnw:322-324
###################################################
(P <- pureboost(L))  # pure boost
P - t(P)   # check for symmetry


###################################################
### code chunk number 9: lorentz.Rnw:329-333
###################################################
(U <- orthog(L))                  # rotation matrix
U[2:4,2:4]                        # inspect the spatial components
round(crossprod(U) - diag(4),10)  # check for orthogonality
## zero to within numerical uncertainty


###################################################
### code chunk number 10: lorentz.Rnw:343-345
###################################################
sol(299792458)
sol()


###################################################
### code chunk number 11: lorentz.Rnw:353-355
###################################################
u <- as.3vel(c(100,200,300))
as.4vel(u)


###################################################
### code chunk number 12: lorentz.Rnw:362-363
###################################################
gam(u)


###################################################
### code chunk number 13: lorentz.Rnw:368-369
###################################################
gam(u)-1


###################################################
### code chunk number 14: lorentz.Rnw:375-376
###################################################
gamm1(u)


###################################################
### code chunk number 15: lorentz.Rnw:381-382
###################################################
boost(u)


###################################################
### code chunk number 16: lorentz.Rnw:396-398
###################################################
sol(299792458)
disp  <- c(1,1,0,0)


###################################################
### code chunk number 17: lorentz.Rnw:406-407
###################################################
ptm(to_natural=TRUE) %*% disp


###################################################
### code chunk number 18: lorentz.Rnw:417-418
###################################################
ptm(to_natural=TRUE,change_time=FALSE) %*% disp


###################################################
### code chunk number 19: lorentz.Rnw:426-429
###################################################
sol(1)
B1 <- boost((2:4)/10) %*% boost(c(-5,1,3)/10)
orthog(B1)[2:4,2:4]


###################################################
### code chunk number 20: lorentz.Rnw:436-439
###################################################
sol(10)
B2 <- boost(2:4) %*% boost(c(-5,1,3))  # exactly the same as B1 above
orthog(B2)[2:4,2:4]


###################################################
### code chunk number 21: setcinfinite
###################################################
sol(Inf)


###################################################
### code chunk number 22: parallelograminfinitec
###################################################
u <- as.3vel(1:3)
v <- as.3vel(c(-6,8,3))
u+v
v+u


###################################################
### code chunk number 23: boostinfc
###################################################
boost(u)
boost(u+v)
boost(u) %*% boost(v)


###################################################
### code chunk number 24: gamm1google
###################################################
gamm1(1e100)


###################################################
### code chunk number 25: rboostinfc
###################################################
set.seed(0)
options(digits=3)
(B <- rboost(1))  # random boost, speed 1


###################################################
### code chunk number 26: orthoginfc
###################################################
orthog(B)
pureboost(B)


###################################################
### code chunk number 27: infinitelightspeedspatialdisplacement
###################################################
(u <- as.3vel(c(10,0,0))) # velocity of 10, parallel to x axis
(B <- boost(u))
d <- c(0,1,0,0) # displacement of distance one, parallel to the x-axis
B %*% d


###################################################
### code chunk number 28: classicalboostspatial
###################################################
d <- c(1,0,0,0) # displacement of one unit of time, no spatial component
B %*% d


###################################################
### code chunk number 29: lorentz.Rnw:533-538
###################################################
sol(10)
u <- as.3vel(c(5,-6,4))
(U <- as.4vel(u))
B <- boost(U)
B %*% as.vector(U)


###################################################
### code chunk number 30: lorentz.Rnw:550-553
###################################################
u <- 1:7  # speed in the x-direction [c=10]
jj <- cbind(gam(u),gam(u)*u,0,0)
(U <- as.4vel(jj))


###################################################
### code chunk number 31: lorentz.Rnw:558-559
###################################################
(B <- boost(as.3vel(c(6,0,0))))  # 60% speed of light


###################################################
### code chunk number 32: lorentz.Rnw:569-570
###################################################
U %*% t(B)


###################################################
### code chunk number 33: lorentz.Rnw:575-576
###################################################
is.consistent.4vel(U %*% t(B))


###################################################
### code chunk number 34: lorentz.Rnw:583-584
###################################################
tcrossprod(U,B)


###################################################
### code chunk number 35: lorentz.Rnw:601-602
###################################################
solve(B)


###################################################
### code chunk number 36: lorentz.Rnw:607-608
###################################################
tcrossprod(U,solve(B))


###################################################
### code chunk number 37: lorentz.Rnw:615-616
###################################################
is.consistent.4vel(tcrossprod(U,solve(B)))


###################################################
### code chunk number 38: lorentz.Rnw:625-629
###################################################
sol(100)
B1 <- boost(r3vel(1)) %*% boost(r3vel(1))
B2 <- boost(r3vel(1)) %*% boost(r3vel(1)) 
(U <- r4vel(5))


###################################################
### code chunk number 39: lorentz.Rnw:635-639
###################################################
U %*% t(B1) %*% t(B2)
U %*% t(B2 %*% B1)    # note order of operations
tcrossprod(U, B2 %*% B1)
U %>% tcrossprod(B2 %*% B1)


###################################################
### code chunk number 40: lorentz.Rnw:649-654
###################################################
U %*% B  # Young Frankenstein: Do Not Use This Brain!
## The above idiom is incorrect.  See
## https://www.youtube.com/watch?v=m7-bMBuVmHo&t=1s
## (in particular @1:08) for a technical explanation of why 
## this is a Very Bad Idea (tm).


###################################################
### code chunk number 41: lorentz.Rnw:669-672
###################################################
sol(1)        # revert to natural units 
D <- dust(1)  # Dust is the simplest nontrivial SET, with 
D             # only one nonzero component


###################################################
### code chunk number 42: lorentz.Rnw:679-681
###################################################
B <- boost(as.3vel(c(0.0,0.8,0.0)))
transform_uu(D,B)


###################################################
### code chunk number 43: lorentz.Rnw:692-695
###################################################
pg <- photongas(3)
pg
transform_uu(pg,B)


###################################################
### code chunk number 44: lorentz.Rnw:704-706
###################################################
raise(transform_dd(lower(pg),lower(B)))
raise(transform_dd(lower(pg),lower(B))) - transform_uu(pg,B) #zero to numerical precision


###################################################
### code chunk number 45: lorentz.Rnw:718-724
###################################################
B1 <- boost(as.3vel(c(0.5,-0.4,0.6)))
B2 <- boost(as.3vel(c(0.1,-0.1,0.3)))
pf <- perfectfluid(4,1)
pf
pf %>% transform_uu(B1) %>% transform_uu(B2)
pf %>% transform_uu(B2 %*% B1)  # should match


###################################################
### code chunk number 46: lorentz.Rnw:730-731
###################################################
lower(pf) %>% transform_dd(lower(B1) %*% lower(B2)) %>% raise()


###################################################
### code chunk number 47: lorentz.Rnw:743-746
###################################################
sol(10)
pf_rest <- perfectfluid(1,4)
pf_rest


###################################################
### code chunk number 48: lorentz.Rnw:754-757
###################################################
u <- as.3vel(3:5)
pf_moving <- perfectfluid(1,4,u)
pf_moving


###################################################
### code chunk number 49: lorentz.Rnw:764-765
###################################################
transform_uu(perfectfluid(1,4,u),boost(u))


###################################################
### code chunk number 50: lorentz.Rnw:775-777
###################################################
sol(1)
(A <- as.photon(as.3vel(cbind(0.9,1:5/40,5:1/40))))


###################################################
### code chunk number 51: lorentz.Rnw:785-786
###################################################
inner4(A)


###################################################
### code chunk number 52: lorentz.Rnw:793-794
###################################################
tcrossprod(A,boost(as.3vel(c(0.7,0,0))))


###################################################
### code chunk number 53: lorentz.Rnw:801-802
###################################################
tcrossprod(A,boost(as.3vel(c(-0.7,0,0))))


###################################################
### code chunk number 54: lorentz.Rnw:807-808
###################################################
tcrossprod(A,solve(boost(as.3vel(c(0.7,0,0)))))


###################################################
### code chunk number 55: lorentz.Rnw:824-826
###################################################
m <- c(1,1,1)
B <- boost(as.3vel(c(0.5,0,0)))


###################################################
### code chunk number 56: lorentz.Rnw:836-838
###################################################
A
(A <- as.4mom(A %*% t(B)))


###################################################
### code chunk number 57: lorentz.Rnw:845-846
###################################################
(A <- reflect(A,m))


###################################################
### code chunk number 58: lorentz.Rnw:855-856
###################################################
(A <- as.4mom(A %*% solve(t(B))))


###################################################
### code chunk number 59: lorentz.Rnw:865-867
###################################################
A <- as.photon(as.3vel(cbind(0.9,1:5/40,5:1/40)))
A %>% tcrossprod(B) %>% reflect(m) %>% tcrossprod(solve(B)) %>% as.4mom


###################################################
### code chunk number 60: disco_ball
###################################################
dfun <- function(n){matrix(rnorm(n*3),ncol=3) %>% sweep(1, sqrt(rowSums(.^2)),`/`)}
(disco <- dfun(7))


###################################################
### code chunk number 61: lorentz.Rnw:884-886
###################################################
p <- as.photon(c(1,0,0))
reflect(p,disco)


###################################################
### code chunk number 62: disco_reflect_percentage
###################################################
table(reflect(p,dfun(1000))[,2]>0) # should be TRUE with probability sqrt(0.5)


###################################################
### code chunk number 63: relativistic_disco
###################################################
B <- boost(as.3vel(c(0.5,0,0)))
p %>% tcrossprod(B) %>% reflect(disco) %>% tcrossprod(solve(B))


###################################################
### code chunk number 64: lorentz.Rnw:924-935
###################################################
sol(1)
light_start <- as.photon(as.3vel(cbind(0.9,1:5/40,5:1/40)))
m <- c(1,0,0)     # mirror normal to x-axis
B1 <- boost(as.3vel(c(-0.5, 0.1, 0.0)))
B2 <- boost(as.3vel(c( 0.2, 0.0, 0.0)))
B3 <- boost(as.3vel(c( 0.0, 0.0, 0.6)))
B <- B1 %*% B2 %*% B3   # matrix multiplication is associative!
light <- light_start %*% t(B)
light <- reflect(light,m)
light <- as.4mom(light %*% solve(t(B)))
light


###################################################
### code chunk number 65: lorentz.Rnw:943-944
###################################################
light_start %>% tcrossprod(B) %>% reflect(m) %>% tcrossprod(solve(B)) %>% as.4mom


###################################################
### code chunk number 66: lorentz.Rnw:948-956
###################################################
sol(Inf)
light_start <- as.photon(as.3vel(cbind(0.9,1:5/40,5:1/40)))
B1 <- boost(as.3vel(c(-0.5, 0.1, 0.0)))
B2 <- boost(as.3vel(c( 0.2, 0.0, 0.0)))
B3 <- boost(as.3vel(c( 0.0, 0.0, 0.6)))
B <- B1 %*% B2 %*% B3
light_start
light_start %>% tcrossprod(B) %>% reflect(m) %>% tcrossprod(solve(B)) %>% as.4mom


###################################################
### code chunk number 67: kickoff
###################################################
sol(1)
u <- as.3vel(c(-0.7,+0.2,-0.3))
v <- as.3vel(c(+0.3,+0.3,+0.4))
w <- as.3vel(c(+0.1,+0.3,+0.8))
x <- as.3vel(c(-0.2,-0.1,-0.9))
u


###################################################
### code chunk number 68: try
###################################################
u+v
v+u


###################################################
### code chunk number 69: lorentz.Rnw:1028-1029
###################################################
(u+v)-gyr(u,v,v+u)


###################################################
### code chunk number 70: funcid
###################################################
f <- gyrfun(u,v)
(u+v)-f(v+u)    # should be zero


###################################################
### code chunk number 71: vec9
###################################################
u9 <- r3vel(9)
u9


###################################################
### code chunk number 72: vecfun
###################################################
f <- gyrfun(u9,v)
f(x)


###################################################
### code chunk number 73: u9+x
###################################################
u9+x


###################################################
### code chunk number 74: nonass
###################################################
(u+v)+w
u+(v+w)


###################################################
### code chunk number 75: nonass1
###################################################
(u+(v+w)) - ((u+v)+gyr(u,v,w))
((u+v)+w) - (u+(v+gyr(v,u,w)))


###################################################
### code chunk number 76: viss
###################################################
u <- as.3vel(c(0.4,0,0))
v <- seq(as.3vel(c(0.4,-0.2,0)), as.3vel(c(-0.3,0.9,0)),len=20)
w <- as.3vel(c(0.8,-0.4,0))


###################################################
### code chunk number 77: comfail1_fig
###################################################
comm_fail1(u=u, v=v)


###################################################
### code chunk number 78: comfail2_fig
###################################################
comm_fail2(u=u, v=v)


###################################################
### code chunk number 79: assfail_fig
###################################################
ass_fail(u=u, v=v, w=w, bold=10)


###################################################
### code chunk number 80: defuvw
###################################################
 u <- as.3vel(c(+0.5,0.1,-0.2))
 v <- as.3vel(c(+0.4,0.3,-0.2))
 w <- as.3vel(c(-0.3,0.2,+0.2))


###################################################
### code chunk number 81: lorentz.Rnw:1162-1165
###################################################
jj1 <- u %>% add(v)
jj2 <- u+v
speed(jj1-jj2)


###################################################
### code chunk number 82: lorentz.Rnw:1170-1173
###################################################
jj1 <- u %>% add(v) %>% add(w)
jj2 <- (u+v)+w
speed(jj1-jj2)


###################################################
### code chunk number 83: lorentz.Rnw:1180-1183
###################################################
jj1 <- u %>% add(v %>% add(w))
jj2 <- u+(v+w)
speed(jj1-jj2)


###################################################
### code chunk number 84: funcnotation
###################################################
x <- as.3vel(c(0.7, 0.0, -0.7))
y <- as.3vel(c(0.1, 0.3, -0.6))
u <- as.3vel(c(0.0, 0.8, +0.1))   # x,y,u: single three-velocities
v <- r3vel(5,0.9)
w <- r3vel(5,0.8)                 # v,w: vector of three-velocities
f <- gyrfun(u,v)
g <- gyrfun(v,u)


###################################################
### code chunk number 85: testeq3to8
###################################################
max(speed((u+v) - f(v+u)))              # equation 3
max(abs(prod3(f(x),f(y)) - prod3(x,y))) # equation 4
max(speed(f(x+y) - (f(x)+f(y))))        # equation 5
max(speed(f(g(x)) - g(f(x))))           # equation 6
max(speed((u+(v+w)) - ((u+v)+f(w))))    # equation 7
max(speed(((u+v)+w) - (u+(v+g(w)))))    # equation 8


