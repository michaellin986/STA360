for(i in 1:n){
  if(xval[i]<0){
    xval[i]=0
  }
  
  if(xval[i]>1){
    xval[i]=1
  }
  
  if(yval[i]<0){
    yval[i]=0
  }
  
  if(yval[i]>1){
    yval[i]=1
  }
}