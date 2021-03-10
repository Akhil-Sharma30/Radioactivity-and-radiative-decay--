Ui=input('Enter intial number of atoms > ');
tau=input('Define decay constant > ');
dt=input('Define time increment > ');

Uarray=zeros(1,100);
t=zeros(1,100)* dt;

Uarray(1)=Ui;

for i=1:99
    Uarray(i+1)=Uarray(i)-dt*(Uarray(i)/tau);
    t(i+1)=t(i)+dt;
end

analytical=Ui*exp(-t/tau);

figure(1)
plot(t,Uarray,'blue'); %Plotting numerical solution 
hold on 
scatter(t,Ui*exp(-t/tau),'red'); %Plotting analytical solution
title('Radioacive Decay Rate Using Euler Method')
xlabel('Time(s)');
ylabel('Number of Nuclei');
legend('Euler Method','Analytical Solution');
hold off

figure(2) %Plotting difference graph
plot(t, Uarray-Ui*exp(-t/tau));
hold on
title('Difference Graph');
xlabel('Time');
ylabel('Error Percentage');
hold off
