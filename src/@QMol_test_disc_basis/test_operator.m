function test_operator(obj)
%test_operator unit tests for the operator overload in the class

    % Initialization
    obj.showSection('Object comparison');
    X                   =   (-15:.1:20);
    V                   =   [exp(-(X-2).^2),exp(-(X-1).^2/.7),exp(-X.^2),exp(-(X+1).^2/2)];
    d_1                 =   QMol_disc_basis('xspan',X,'basis',V);
    d_2                 =   QMol_disc_basis('xspan',X + .5e-10,'basis',V+.91e-10);

    % ==
    obj.showResult('== (eq operator)',d_1 == d_2);

    % ~=
    d_2.set('xspan',X + 1.01e-10);
    obj.showResult('~= (value mismatch)',d_1 ~= d_2);

    % ~=
    d_2.set('xspan',[]);
    obj.showResult('~= (size mismatch)',d_1 ~= d_2);
end