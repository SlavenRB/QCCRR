
qccr(x,y,inc)

qccr <- function(x, y, inc) {

#--------------------------------------------------------------------------.
#	1`	dimenzije ulaznih matrica; inter i kros korelacije
#--------------------------------------------------------------------------.
nx <- ncol(x)
ny <- ncol(y)
cx <- cor(x)
cy <- cor(y)
cxy <- cor(x,y)
ns <- nrow(x)

#--------------------------------------------------------------------------.
#	2`	dekompozicija Cholesk-ok; racunanje inverza trougaonih matrica
#--------------------------------------------------------------------------.
r1 <-chol(cx)
r2 <-chol(cy)
r1inv <- solve(r1)
r2inv <- solve(r2)

#--------------------------------------------------------------------------.
#	3 izracunavanje omega matrice i njena singularna dekompozicija
#--------------------------------------------------------------------------.
if(nx <= ny) {
	omega <- t(r1inv) %*% cxy  %*% r2inv
} else {
	omega <- t(r2inv) %*% t(cxy)  %*% r1inv
}
domega <-svd(omega) 
dlam <- matrix (domega$d , ncol=1)

#--------------------------------------------------------------------------.
#	4 Kanonicke korelacije i testiranje njihove znacajnosti
#--------------------------------------------------------------------------.
eign <- 1 / (1 - dlam * dlam) - 1
wlam <- 1 / (1+ eign)
n <- nrow(wlam)
wilk <- wlam
df <- wlam
sig <- wlam
bart2 <- wlam
tem <- 1

for(i in 1:n){
	tem <- tem * wlam [n - i + 1]
	df[i] <- (nx - n + i) * (ny - n + i)
	bart2[i] <-  - (ns - (nx + ny + 3)/ 2) * logb(tem) 
	wilk[i] <- tem
}

bart2 <- bart2[order(-bart2),] 
bart2 <- matrix(bart2, ncol=1)
df <- df[order(-df)] 
df <- matrix(df, ncol=1)

for(i in 1:n){
	chi <- bart2 [n - i + 1]
	dof <- df[n - i + 1]
	sig[i] <- 1 - pchisq(chi, dof)
}

dlam <- dlam[order(dlam)] 
dlam <- matrix(dlam, ncol=1)
spoj <- cbind(dlam, sig, wilk)
spoj <- spoj[order(-dlam),]
test <- cbind ((spoj[,1]), (spoj[,3]), bart2, df, (spoj[ ,2]))
test <- round(test, digits=3)
colnames(test) <- c("Rho", "Lambda", "Hi2", "df", "sig")

#-------------------------------------------------------------------------.
#	5 Zadrzavanje znacajnih kanonickih korelacija
#-------------------------------------------------------------------------.
znac <-0
for (i in sig){
	if (i < inc){
	znac <- znac + 1
	}
}

#---------------------------------------------------------------------------.
#	6 Izracunavanje kanonickih koeficijenata levog skupa varijabli
#---------------------------------------------------------------------------.

if (znac > 0){
	if(nx <= ny) {
		a <- r1inv %*% domega$u
	} else {
		a <- r1inv %*% domega$v
	}
	if(ny < nx){
		a <- a[ ,1:ny]
	}
	aznac <- a[ ,1:znac]

#---------------------------------------------------------------------------.
#	7 Izracunavanje kanonickih koeficijenata desnog skupa varijabli         
#---------------------------------------------------------------------------.
	if(nx <= ny) {
		b <- r2inv %*% domega$v
	} else {
		b <- r2inv %*% domega$u
	}
	if(nx <= ny){
		b <- b[ ,1:nx]
	} 
	bznac <- b[ ,1:znac]

#---------------------------------------------------------------------------.
#	8 Izracunavanje struktura kanonickih faktora levog skupa varijabli
#---------------------------------------------------------------------------.
	tem1 <- cx %*% a
	tem1z <- tem1[ ,1:znac]
	tem1z <- matrix(tem1z, ncol=znac)
	rownames(tem1z)<- names(x)
	tem1z <- round(tem1z, digits=3)

#---------------------------------------------------------------------------.
#	9 Izracunavanje struktura kanonickih faktora desnog skupa varijabli
#---------------------------------------------------------------------------.
	tem2 <- cy %*% b
	tem2z <- tem2[ ,1:znac]
	tem2z <- matrix(tem2z, ncol=znac)
	rownames(tem2z)<- names(y)
	tem2z <- round(tem2z, digits=3)

#---------------------------------------------------------------------------.
#	10 Izracunavanje kanonickih kros-faktora levog skupa varijabli
#---------------------------------------------------------------------------.
	tem12 <- cxy %*% b
	tem12z <- tem12[ ,1:znac]
	tem12z <- matrix(tem12z, ncol=znac)
	rownames(tem12z)<- names(x)
	tem12z <- round(tem12z, digits=3)

#---------------------------------------------------------------------------.
#	11 Izracunavanje kanonickih kros-faktora levog skupa varijabli
#---------------------------------------------------------------------------.
	tem21 <- t(cxy) %*% a
	tem21z <- tem21 [ ,1:znac]
	tem21z <- matrix(tem21z, ncol=znac)
	rownames(tem21z)<- names(y)
	tem21z <- round(tem21z, digits=3)

#----------------------------------------------------------------------------.      
# 	12 Relativne varijanse, indeksi prepokrivanja i mere generalizabilnosti levog skupa varijabli                                
#----------------------------------------------------------------------------.  
	cst1 <- colSums(tem1z * tem1z)
	cst1x <- cst1 /nx
	jedanc <- matrix (rep(1, znac), ncol=1)
	alpha1 <- (nx/(nx-1)) * (jedanc - jedanc / cst1)
	cs3 <- (colSums(tem12z * tem12z))/nx
	cs3 <- matrix(cs3, ncol=1)
	test1 <- cbind (cst1x, cs3, alpha1)  
	test1 <- round(test1, digits=3)
	colnames(test1) <- c("Var.", "Prepok.", "General.")

#----------------------------------------------------------------------------.      
# 	13 Relativne varijanse, indeksi prepokrivanja i mere generalizabilnosti desnog skupa varijabli                                
#----------------------------------------------------------------------------.
	cst2 <- colSums(tem2z * tem2z)
	cst2y <- cst2 /ny
	alpha2 <- (ny/(ny-1)) * (jedanc - jedanc / cst2)
	cs4 <- (colSums(tem21z * tem21z))/ny
	cs4 <- matrix(cs4, ncol=1)
	test2 <- cbind (cst2y, cs4, alpha2)  
	test2 <- round(test2, digits=3)
	colnames(test2) <- c("Var.", "Prepok.", "General.")

	aznac <- matrix(aznac, ncol=1)
	bznac <- matrix(bznac, ncol=1)
	aznac <- round(aznac, digits=3)
	bznac <- round(bznac, digits=3)
}

# KANONICKA ANALIZA KOVARIJANSI

#---------------------------------------------------------------------------. 
#	14 Singularna dekompozicija matrice kroskorelacija                          
#---------------------------------------------------------------------------. 
dxy <- svd(cxy)

#---------------------------------------------------------------------------. 
#	15 Zadrzavaju se singularne vrednosti matrice r12 koje su vece od 
#	aritmeticke sredine svih njenih singularnih vrednosti                          
#---------------------------------------------------------------------------. 
as <- mean(dxy$d)
t <-(dxy$d >= as)
znacq <- length(dxy$d[t])
txy <- svd((cor(x,y)),nu=znacq, nv=znacq)
x1z <- dxy$u[ ,1:znacq]
x2z <- dxy$v[ ,1:znacq]

#---------------------------------------------------------------------------. 
#	16 Izracunavanje kvazikanonickih koeficijenata, sklopa, strukture,    
#	interkorelacija zadrzanih kvazikanonickih faktora 
#	i testiranje njihove znacajnosti                           
#---------------------------------------------------------------------------. 
c <- t(txy$u) %*% cxy %*% txy$v
mm1 <- (t(txy$u)) %*% cx %*% txy$u
mm2 <- (t(txy$v)) %*% cy %*% txy$v

x1z <- round (x1z , digits=3)
x2z <- round (x2z , digits=3)


if (znacq > 1){
	d1 <- diag(mm1)
	d1 <- sqrt(d1)
	d1 <- 1/d1
	d1 <- diag(d1)
	d2 <- diag(mm2)
	d2 <- sqrt(d2)
	d2 <- 1/d2
	d2 <- diag(d2)
	m1 <- d1 %*% mm1 %*% d1
	m2 <- d2 %*% mm2 %*% d2
	invm1 <- solve(m1)
	invm2 <- solve(m2)
	f1 <- cx %*% txy$u %*% d1
	f2 <- cy %*% txy$v %*% d2
	a1 <- f1 %*% invm1
	a2 <- f2 %*% invm2
	ro <- d1 %*% c %*% d2
	ro <- diag(ro)
	ro2 <- ro * ro
	f12 <- cxy %*% txy$v %*% d2
	f21 <- t(cxy) %*% txy$u %*% d1
	df1 <- 1
	df2 <- ns - 2
	ftest1 <- matrix (rep(1, znacq), ncol=1)
	ftest2 <- ftest1
	sig1 <- ftest1
	for(i in 1:znacq){
		ftest1[i] <- ro2[i] * ((ns - 2)/(1 - ro2[i])) 
		ftest2 <-ftest1[i]
		sig1[i] <- 1 - pf(ftest2, df1, df2)
	}
	ro <- round(ro, digits=3)
	ro2 <- round(ro2, digits=3)
	ftest1 <- round(ftest1, digits=3)
	sig1 <- round(sig1, digits=5)
	lav <- cbind(ro, ro2, ftest1,sig1)
	colnames(lav) <- c("ro", "ro2", "f-test", "sig")

#---------------------------------------------------------------------------.             
#	17 Relativne varijanse, indeksi prepokrivanja i mere generalizabilnosti                  
#	levog skupa varijabli                                                                 
#---------------------------------------------------------------------------.                     
	jedan <- matrix (rep(1, znacq), ncol=1)
	csf1 <- colSums(f1 * f1)
	red1 <- csf1 /nx
	red1 <- matrix(red1, ncol=1)
	red3 <- colSums(f12 * f12)/nx
	red3 <- matrix(red3, ncol=1)
	csf1 <- matrix(csf1, ncol=1)
	beta1 <- (nx /(nx - 1)) * (jedan - jedan / csf1)  
	test1c <- cbind (red1, red3, beta1)  
	test1c <- round(test1c, digits=3)
	colnames(test1c) <- c("Var.", "Prepok.", "General.")

#---------------------------------------------------------------------------.             
#	18 Relativne varijanse, indeksi prepokrivanja i mere generalizabilnosti                  
# 	levog skupa varijabli                                                                 
#---------------------------------------------------------------------------. 
	csf2 <- colSums(f2 * f2)
	red2 <- csf2 /ny
	red2 <- matrix(red2, ncol=1)
	red4 <- colSums(f21 * f21)/ny
	red4 <- matrix(red4, ncol=1)
	csf2 <- matrix(csf2, ncol=1)
	beta2 <- (ny /(ny - 1)) * (jedan - jedan / csf2)  
	test2c <- cbind (red2, red4, beta2)  
	test2c <- round(test2c, digits=3)
	colnames(test2c) <- c("Var.", "Prepok.", "General.")	
	a1 <- round (a1 , digits=3)
	m1 <- round (m1 , digits=3)
	f1 <- round (f1 , digits=3)
	f12 <- round (f12 , digits=3)
	a2 <- round (a2 , digits=3)
	m2 <- round (m2 , digits=3)
	f2 <- round (f2 , digits=3)
	f21 <- round (f21 , digits=3)
	
#---------------------------------------------------------------------------.             
#	19 U slucaju da je izdvojena samo jedna znacajna kvazikanonikcka funkcija                                                                                   
#---------------------------------------------------------------------------.
	
}	else if (znacq == 1){

	d1 <- sqrt(mm1)
	d2 <- sqrt(mm2)
	d1 <- 1/d1
	d2 <- 1/d2	
	f1 <- cx %*% txy$u %*% d1
	f2 <- cy %*% txy$v %*% d2
	f12 <- cxy %*% txy$v %*% d2
	f21 <- t(cxy) %*% txy$u %*% d1	
	ro <- d1 %*% c %*% d2
	ro2 <- ro * ro
	df1 <- 1
	df2 <- ns - 2
	ftest <- ro2*(df2/(1-ro2))
	sig1q <- 1 - pf(ftest, 1, df2) 
	ro <- round(ro, digits=3)
	ro2 <- round(ro2, digits=3)
	ftest <- round(ftest, digits=3)
	sig1q <- round(sig1q, digits=5)
	lav <- cbind(ro, ro2, ftest, df1, df2, sig1q)
	colnames(lav) <- c("ro", "ro2", "f-test", "df1", "df2", "sig")
		if (sig1q < inc){
		desnilav <- matrix (c(txy$u, f1, f12), ncol=3)
		colnames(desnilav) <- c("X","F","F1")
		rownames(desnilav)<- names(x)
		desnilav <- round(desnilav, digits=3)
		levilav <- matrix (c(txy$v, f2, f21), ncol=3)
		colnames(levilav) <- c("X","F","F2")
		rownames(levilav)<- names(y)
		levilav <- round(levilav, digits=3)
		}	
	
	jedan <- matrix (rep(1, znacq), ncol=1)
	csf1 <- colSums(f1 * f1)
	red1 <- csf1 /nx
	red1 <- matrix(red1, ncol=1)
	red3 <- colSums(f12 * f12)/nx
	red3 <- matrix(red3, ncol=1)
	csf1 <- matrix(csf1, ncol=1)
	beta1 <- (nx /(nx - 1)) * (jedan - jedan / csf1)  
	test1c <- cbind (red1, red3, beta1)  
	test1c <- round(test1c, digits=3)
	colnames(test1c) <- c("Var.", "Prepok.", "General.")
		
	csf2 <- colSums(f2 * f2)
	red2 <- csf2 /ny
	red2 <- matrix(red2, ncol=1)
	red4 <- colSums(f21 * f21)/ny
	red4 <- matrix(red4, ncol=1)
	csf2 <- matrix(csf2, ncol=1)
	beta2 <- (ny /(ny - 1)) * (jedan - jedan / csf2)  
	test2c <- cbind (red2, red4, beta2)  
	test2c <- round(test2c, digits=3)
	colnames(test2c) <- c("Var.", "Prepok.", "General.")		
}


if (znac > 0){

#---------------------------------------------------------------------------.            
#	20 Izracunavanje korelacija izmedju kanonickih i kvazikanonickih faktora 
#	u levom skupu varijabli                                                    
#---------------------------------------------------------------------------.     
	aznac <-matrix(aznac, ncol=znac)
	cq1 <- t(aznac) %*% f1
	cq1 <- round(cq1, digits=4)

#---------------------------------------------------------------------------.            
#	21 Izracunavanje korelacija izmedju kanonickih i kvazikanonickih faktora 
#	u desnom skupu varijabli                                                    
#---------------------------------------------------------------------------.     
	bznac <-matrix(bznac, ncol=znac)
	cq2 <- t(bznac) %*% f2
	cq2 <- round(cq2, digits=4)

#-------------------------------------------------------------------------.
#	22 Izracunavanje i Tucker-ovih koeficijenata kongruencije
#	kanonickih i kvazikanonickih faktora
#--------------------------------------------------------------------------.

	g1 <- t(f1) %*% f1
	g2 <- t(f2) %*% f2
		if(znacq > 1){
			g1 <- diag(g1)
			g1 <- sqrt(g1)
			g1 <- 1/g1
			g1 <- diag(g1)
			g2 <- diag(g2)
			g2 <- sqrt(g2)
			g2 <- 1/g2
			g2 <- diag(g2)
		} else {
			g1 <- sqrt(g1)
			g1 <- 1/g1
			g2 <- sqrt(g2)
			g2 <- 1/g2
		}
	h1 <- t(tem1z) %*% tem1z
	h2 <- t(tem2z) %*% tem2z
	
		if (znac >1){
			h1 <- diag(h1)
			h1 <- sqrt(h1)
			h1 <- 1/h1
			h1 <- diag(h1)
			h2 <- diag(h2)
			h2 <- sqrt(h2)
			h2 <- 1/h2
			h2 <- diag(h2)
		} else {
			h1 <- sqrt(h1)
			h1 <- 1/h1
			h2 <- sqrt(h2)
			h2 <- 1/h2
		}
	kong1 <- t(tem1z) %*% f1
	kong1 <- h1 %*% kong1 %*% g1
	kong2 <- t(tem2z) %*% f2
	kong2 <- h2 %*% kong2 %*% g2
	
	kong1 <- round (kong1 , digits=4)
	kong2 <- round (kong2 , digits=4)
}

cxr <- round(cx, digits=3)
cyr <- round(cy, digits=3)
cxyr <- round(cxy, digits=3)

#----------------------------------------------------------------------------.      
#	23 Stampanje rezultata                            
#----------------------------------------------------------------------------.
if(znac > 0 & znacq > 1){
	list (
	"Interkorelacije levog skupa varijabli" = cxr,
	"Interkorelacije desnog skupa varijabli" = cyr,
	"Kroskorelacije levog i desnog skupa varijabli" = cxyr,
	"Koeficijenti kanonickih korelacija i njihova znacajnost:" = test, 
	"Kanonicki koeficijenti levog skupa varijabli" = aznac,
	"Kanonicki koeficijenti desnog skupa varijabli" = bznac,
	"Kanonicki faktori levog skupa varijabli" = tem1z,
	"Kanonicki faktori desnog skupa varijabli" = tem2z,
	"Kanonicki kros-faktori levog skupa varijabli" = tem12z, 
	"Kanonicki kros-faktori desnog skupa varijabli" = tem21z,
	"Analiza prepokrivanja levog skupa"=test1, 
	"Analiza prepokrivanja desnog skupa"=test2,
	"Kvazikanonicke korelacije i testovi znacajnosti" = lav,
	"Kvazikanonicki koeficijenti levog skupa varijabli" = x1z,
	"Sklop levog skupa varijabli" = a1,
	"Interkorelacije levog skupa kvazikanonickih faktora" = m1,
	"Struktura levog skupa varijabli kvazikanonickih faktora"= f1,
	"Kros-struktura levog skupa" =f12,
	"Kvazikanonicki koeficijenti desnog skupa varijabli" = x2z,
	"Sklop levog skupa varijabli" = a2,
	"Interkorelacije levog skupa kvazikanonickih faktora" = m2,
	"Struktura levog skupa varijabli kvazikanonickih faktora"= f2,
	"Kros-struktura levog skupa" =f21,
	"Analiza prepokrivanja levog skupa" = test1c,
	"Analiza prepokrivanja desnog skupa" = test2c,
	"Korelacije levog skupa kanonickih i kvazikanonickih faktora" = cq1,
	"Korelacije desnog skupa kanonickih i kvazikanonickih faktora" = cq2,
	"Kongruencija levog skupa kanonickih i kvazikanonickih faktora" = kong1,
	"Kongruencija desnog skupa kanonickih i kvazikanonickih faktora" = kong2
	)
	} else if (znac > 0 & znacq == 1 & sig1q < inc){
	list (
	"Interkorelacije levog skupa varijabli" = cxr,
	"Interkorelacije desnog skupa varijabli" = cyr,
	"Kroskorelacije levog i desnog skupa varijabli" = cxyr,
	"Koeficijenti kanonickih korelacija i njihova znacajnost:" = test, 
	"Kanonicki koeficijenti levog skupa varijabli" = aznac,
	"Kanonicki koeficijenti desnog skupa varijabli" = bznac,
	"Kanonicki faktori levog skupa varijabli" = tem1z,
	"Kanonicki faktori desnog skupa varijabli" = tem2z,
	"Kanonicki kros-faktori levog skupa varijabli" = tem12z, 
	"Kanonicki kros-faktori desnog skupa varijabli" = tem21z,
	"Analiza prepokrivanja levog skupa"=test1, 
	"Analiza prepokrivanja desnog skupa"=test2,
	"Kvazikanonicke korelacije i testovi znacajnosti" = lav,
	"Koeficijenti, struktura i krosstruktura levog skupa varijabli" = levilav,
	"Koeficijenti, struktura i krosstruktura desnog skupa varijabli" = desnilav,
	"Analiza prepokrivanja levog skupa" = test1c, 
	"Analiza prepokrivanja desnog skupa" = test2c	
	)
	} else if (znac > 0 & znacq == 1 & sig1q >= inc){
	list (
	"Interkorelacije levog skupa varijabli" = cxr,
	"Interkorelacije desnog skupa varijabli" = cyr,
	"Kroskorelacije levog i desnog skupa varijabli" = cxyr,
	"Koeficijenti kanonickih korelacija i njihova znacajnost:" = test, 
	"Kanonicki koeficijenti levog skupa varijabli" = aznac,
	"Kanonicki koeficijenti desnog skupa varijabli" = bznac,
	"Kanonicki faktori levog skupa varijabli" = tem1z,
	"Kanonicki faktori desnog skupa varijabli" = tem2z,
	"Kanonicki kros-faktori levog skupa varijabli" = tem12z, 
	"Kanonicki kros-faktori desnog skupa varijabli" = tem21z,
	"Analiza prepokrivanja levog skupa"=test1, 
	"Analiza prepokrivanja desnog skupa"=test2,
	"Kvazikanonicke korelacije i testovi znacajnosti" = lav,
	"Analiza prepokrivanja levog skupa" = test1c, 
	"Analiza prepokrivanja desnog skupa" = test2c	
	)
	} else if (znac <= 0 & znacq > 1){
	list ("Interkorelacije levog skupa varijabli" = cxr,
	"Interkorelacije desnog skupa varijabli" = cyr,
	"Kroskorelacije levog i desnog skupa varijabli" = cxyr,
	"Kvazikanonicke korelacije i testovi znacajnosti" = lav,
	"Kvazikanonicki koeficijenti levog skupa varijabli" = x1z,
	"Sklop levog skupa varijabli" = a1,
	"Interkorelacije levog skupa kvazikanonickih faktora" = m1,
	"Struktura levog skupa varijabli kvazikanonickih faktora"= f1,
	"Kros-struktura levog skupa" =f12,
	"Kvazikanonicki koeficijenti desnog skupa varijabli" = x2z,
	"Sklop levog skupa varijabli" = a2,
	"Interkorelacije levog skupa kvazikanonickih faktora" = m2,
	"Struktura levog skupa varijabli kvazikanonickih faktora"= f2,
	"Kros-struktura levog skupa" =f21,
	"Analiza prepokrivanja levog skupa" = test1c,
	"Analiza prepokrivanja desnog skupa" = test2c
	)		
	} else if (znac <= 0 & znacq == 1 & sig1q < inc){
	list ("Interkorelacije levog skupa varijabli" = cxr,
	"Interkorelacije desnog skupa varijabli" = cyr,
	"Kroskorelacije levog i desnog skupa varijabli" = cxyr,
	"Kvazikanonicke korelacije i testovi znacajnosti" = lav,
	"Koeficijenti, struktura i krosstruktura levog skupa varijabli" = levilav,
	"Koeficijenti, struktura i krosstruktura desnog skupa varijabli" = desnilav,
	"Analiza prepokrivanja levog skupa" = test1c, 
	"Analiza prepokrivanja desnog skupa" = test2c	
	)		
	} else if (znac <= 0 & znacq == 1 & sig1q >= inc){
	list ("Interkorelacije levog skupa varijabli" = cxr,
	"Interkorelacije desnog skupa varijabli" = cyr,
	"Kroskorelacije levog i desnog skupa varijabli" = cxyr,
	"Kvazikanonicke korelacije i testovi znacajnosti" = lav,
	"Analiza prepokrivanja levog skupa" = test1c, 
	"Analiza prepokrivanja desnog skupa" = test2c	
	)	
	} else if (znac > 0 & znacq <= 0){
	list (
	"Interkorelacije levog skupa varijabli" = cxr,
	"Interkorelacije desnog skupa varijabli" = cyr,
	"Kroskorelacije levog i desnog skupa varijabli" = cxyr,
	"Koeficijenti kanonickih korelacija i njihova znacajnost:" = test, 
	"Kanonicki koeficijenti levog skupa varijabli" = aznac,
	"Kanonicki koeficijenti desnog skupa varijabli" = bznac,
	"Kanonicki faktori levog skupa varijabli" = tem1z,
	"Kanonicki faktori desnog skupa varijabli" = tem2z,
	"Kanonicki kros-faktori levog skupa varijabli" = tem12z, 
	"Kanonicki kros-faktori desnog skupa varijabli" = tem21z,
	"Analiza prepokrivanja levog skupa"=test1, 
	"Analiza prepokrivanja desnog skupa"=test2
	)
	} 
}



