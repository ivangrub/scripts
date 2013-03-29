%Pull out the genes that are primarily associated with each of the
%Principal Components

%The determining factors will be whether or not the principal component
%scores are z >= 3 away from the origin and have a correlation coefficient
%of at least 0.8 with one of the components.
A = zeros(length(knownGene),3);
B = zeros(length(knownGene),3);

%Parallel to Principal Component Coefficients
for i = 1:length(knownGene)
    A(i,1) = dot(PC(i,:)',coeff(:,1))/(norm(coeff(:,1))*norm(PC(i,:)));
    A(i,2) = dot(PC(i,:)',coeff(:,2))/(norm(coeff(:,2))*norm(PC(i,:)));
    A(i,3) = dot(PC(i,:)',coeff(:,3))/(norm(coeff(:,3))*norm(PC(i,:)));
end

I = (abs(A(:,1)) >= 0.8 & abs(A(:,1)) <= 0.2 & abs(A(:,3)) <= 0.2);

Genes = knownGene(I,2:3);

%Parallel to Principal Component Axis
for i = 1:length(knownGene)
    B(i,1) = dot(PC(i,:),[1 0 0])/(norm(PC(i,:)));
    B(i,2) = dot(PC(i,:),[0 1 0])/(norm(PC(i,:)));
    B(i,3) = dot(PC(i,:),[0 0 1])/(norm(PC(i,:)));
end

I = (abs(B(:,3)) >= 0.8 & abs(B(:,1)) <= 0.2 & abs(B(:,2)) <= 0.2);

Genes2 = knownGene(I,2:3);
