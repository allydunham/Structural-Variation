# Copyright (c) 2017  Genome  Research  Ltd.
# Author: Alistair Dunham
# This  program  is free  software: you  can  redistribute  it and/or  modify  it  under
# the  terms  of the  GNU  General  Public  License  as  published  by the  Free  Software
# Foundation; either  version 3 of the  License , or (at your  option) any  later
# version.
# This  program  is  distributed  in the  hope  that it will be useful , but  WITHOUT
# ANY  WARRANTY; without  even  the  implied  warranty  of  MERCHANTABILITY  or  FITNESS
# FOR A PARTICULAR  PURPOSE. See  the  GNU  General  Public  License  for  more
# details.
# You  should  have  received a copy of the  GNU  General  Public  License  along  with
# this  program. If not , see <http :// www.gnu.org/licenses/>.

## Simple function to create a colour bar in a specified region of a plot

colourBar <- function(fig=c(0.8,1,0,1),colours=c("red","blue"),labels=c(0,0.5,1),main="",lim=c(0,1),srt=90){
  par.backup <- par(no.readonly = TRUE)
  par(fig=fig,new=TRUE,xpd=NA,mar=c(1,1,1,1))
  cols <- colorRampPalette(colours)(100)
  plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab="",ylab="",bty='n')
  sapply(1:100,function(x){rect(0,(x-1)/100,1,x/100,col = cols[x],border = cols[x])})
  axis(4,at = seq(0,1,length.out = length(labels)),labels = labels)
  text(2.5,0.5,labels = main,srt=srt)
  par(par.backup)
}
