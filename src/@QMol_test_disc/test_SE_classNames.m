function test_SE_classNames(obj)
%test_SE_classNames test class names for Schrodinger-equation simulations

    % Component names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('Schrodinger-equation component names');

    obj.showResult('SE_classWaveFunction',      strcmpi(QMol_disc.SE_classWaveFunction,     'QMol_SE_wfcn'));
    obj.showResult('SE_classPotential',         strcmpi(QMol_disc.SE_classPotential,        'QMol_SE_V'));

end