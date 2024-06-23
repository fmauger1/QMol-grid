function test_operator(obj)
%test_operator unit tests for the operator overload in the class

    % Initialization
    obj.showSection('Object comparison');
    X                   =   -15:.1:20;
    d_1                 =   QMol_disc('xspan',X);
    d_2                 =   QMol_disc('xspan',X + .5e-10);

    % ==
    obj.showResult('== (eq operator)',d_1 == d_2);

    % ~=
    d_2.set('xspan',X + 1.01e-10);
    obj.showResult('~= (value mismatch)',d_1 ~= d_2);

    % ~=
    d_2.set('xspan',[]);
    obj.showResult('~= (size mismatch)',d_1 ~= d_2);
end