%Constant parameters
l1 = 1;
l2 = 2;
xe = 0.25;
ye = 2;
x1 = 0;
y1 = 0;
x3 = -1;
y3 = 0;

%Numerical objective function to calculate error of guessed joint values
objectiveFunction = @(angles) [
  xe - (x1 + l1*cos(angles(3)) + l2*cos(angles(3)+angles(1)));
  ye - (y1 + l1*sin(angles(3)) + l2*sin(angles(3)+angles(1)));
  xe - (x3 + l1*cos(angles(4)) + l2*cos(angles(4)+angles(2)));
  ye - (y3 + l1*sin(angles(4)) + l2*sin(angles(4)+angles(2)))
];

%Call a rountine to iteratively update the guess, guiding error to zero
initial_guess = [0;0;0;0];
numerically_solved_angles = fsolve(objectiveFunction, initial_guess);

%Plot the robot
alpha = numerically_solved_angles(3)
beta = numerically_solved_angles(4)

x2 = l1*cos(alpha);
y2 = l1*sin(alpha);

x4 = x3 + l1*cos(beta);
y4 = y3 + l1*sin(beta);

plot([x1;x2;xe;x4;x3],[y1;y2;ye;y4;y3])
hold on
plot([x1;x2;xe;x4;x3],[y1;y2;ye;y4;y3],'or')
hold off
xlabel('x')
ylabel('y')
title('5-bar Closed-chain Robot')
grid on
daspect([1 1 1])