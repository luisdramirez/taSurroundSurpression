%[t1Valid t1Invalid; t2Valid t2Invalid]
Y1 = [avgCollDiffPerTarg(2,1) avgCollDiffPerTarg(3,1); avgCollDiffPerTarg(2,2) avgCollDiffPerTarg(3,2)]; 
bar(Y1)
title('collinear')
legend('valid','invalid')
xlabel('target')
ylabel('average difference contrast')
axis square
ylim([0 1])

Y2 = [avgOrthDiffPerTarg(2,1) avgOrthDiffPerTarg(3,1); avgOrthDiffPerTarg(2,2) avgOrthDiffPerTarg(3,2)];
figure
bar(Y2)
title('orthogonal')
legend('valid','invalid')
xlabel('target')
ylabel('average difference contrast')
axis square
ylim([0 1])


Y3 = [avgBaseDiffPerTarg(2,1) avgBaseDiffPerTarg(3,1); avgBaseDiffPerTarg(2,2) avgBaseDiffPerTarg(3,2)];
figure
bar(Y3)
title('baseline')
legend('valid','invalid')
xlabel('target')
ylabel('average difference contrast')
axis square
ylim([0 1])