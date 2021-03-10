%initializing constants 

q=Q(t);



function f= odeProject(t,y)
f=r*q-{r/{d-r}*t+v}q;
end
 
