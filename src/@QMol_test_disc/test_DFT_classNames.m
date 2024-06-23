function test_DFT_classNames(obj)
%test_DFT_classNames test class names for DFT simulations

    % Component names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('DFT component names');

    obj.showResult('DFT_classDensity',          strcmpi(QMol_disc.DFT_classDensity,          'QMol_DFT_density'));
    obj.showResult('DFT_classOrbital',          strcmpi(QMol_disc.DFT_classOrbital,          'QMol_DFT_orbital'));
    obj.showResult('DFT_classPotential',        strcmpi(QMol_disc.DFT_classPotential,        'QMol_DFT_Vks'));
    obj.showResult('DFT_classPotentialGradient',strcmpi(QMol_disc.DFT_classPotentialGradient,'QMol_DFT_Vks_grad'));

end