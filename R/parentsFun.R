parentsFun <-
function(size, numb){
	man = sample(1:size, numb, replace=TRUE)
	woman = sample(1:size, numb, replace=TRUE)
	cbind(man, woman)
}
