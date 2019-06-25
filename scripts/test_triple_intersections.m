% prototyping a non-svd solution

pLinesind = [.879009  .476806 .0386237
-.707246 .706968 -.0278674
.383775  .923427 .0189677];

l = pLinesind';

l1 = l(:,1);
l2 = l(:,2);
l3 = l(:,3);

l1l1 = l1'*l1;
l1l2 = l1'*l2;
l2l1 = l1l2;
l2l2 = l2'*l2;
l3l1 = l3'*l1;
l3l2 = l3'*l2;

a = cross([l1l1 l2l1 l3l1], [l1l2 l2l2 l3l2]);
