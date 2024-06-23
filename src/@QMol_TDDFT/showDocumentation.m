function showDocumentation(obj)
%showDocumentation displays the documentation reflecting member property 
%   values
    
    % Header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ~obj.isRun
        QMol_doc.showHeader;
    end

    % Time propagation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QMol_doc.showSection('Time propagation');
    fprintf('  * The TDDFT model dynamics matches the canonical Hamiltonian formalism\n');
    fprintf('    discussed in [Mauger 2024].\n')

    ref                 =   ['Mauger 2024',obj.showDoc];

    fprintf('  * Parameters\n');
    fprintf('    Time interval = %s to %s a.u. (%s to %s fs)\n',...
        num2str(obj.T(1),'%5.2f'),                              ...
        num2str(obj.T(end),'%5.2f'),                            ...
        num2str(convertUnit.au2second(obj.T(1))*1e15,'%5.2f'),  ...
        num2str(convertUnit.au2second(obj.T(end))*1e15,'%5.2f'));
    fprintf('    Time step     = %s a.u. (%s as)\n',          ...
        num2str(obj.dt,'%5.3f'),                                ...
        num2str(convertUnit.au2second(obj.dt)*1e18,'%5.3f'));

    if ~isempty(obj.EF),    ref =   [ref,obj.EF.showDocumentation];         end

    if ~isempty(obj.ABC),   ref =   [ref,obj.ABC.showDocumentation];        end

    if obj.sRest,       fprintf('  * Restart option is enabled\n');
    else,               fprintf('  * Restart option is disabled\n');        end

    fprintf('\n');
    
    % Output results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QMol_doc.showSection('Output results');

    if obj.sEDFT   ||   obj.sEKSO
        fprintf('  * DFT and orbital energies:\n'); 
        if obj.sEDFT,   fprintf('    > Save the DFT energy\n');                             end
        if obj.sEKSO,   fprintf('    > Save the Kohn-Sham orbitals'' energy\n');            end
    end

    if obj.sDip   ||   obj.sVel   ||   obj.sAcc
        fprintf('  * Dipole, dipole velocity, and dipole acceleration:\n'); 
        if obj.sDip,    fprintf('    > Save the dipole signal\n');                          if ~isempty(obj.sDipI)
                        fprintf('      + orbital-resolved dipole signal\n');                end, end
        if obj.sVel,    fprintf('    > Save the dipole-velocity signal\n');                 if ~isempty(obj.sVelI)
                        fprintf('      + orbital-resolved dipole-velocity signal\n');       end, end
        if obj.sAcc,    fprintf('    > Save the dipole-acceleration signal\n');             if ~isempty(obj.sAccI)
                        fprintf('      + orbital-resolved dipole-acceleration signal\n');   end, end
    end

    if obj.sIon
        fprintf('  * Ionization statistics:\n    > Save the ionization signal\n');          if ~isempty(obj.sIKSOI)
                        fprintf('      + orbital-resolved ionization signal\n');            end
    end

    if obj.sRho   ||   obj.sKSO   ||   obj.sKSOP
        fprintf('  * Kohn-Sham orbitals and one-body density:\n'); 
        if obj.sRho,    fprintf('    > Save the one-body density\n');                       end
        if obj.sKSO,    fprintf('    > Save the Kohn-Sham orbitals\n');                     end
        if obj.sKSOP,   fprintf('    > Save the projection of the Kohn-Sham orbitals\n');   end
    end

    if isa(obj.sFRho,'function_handle')   ||   isa(obj.sFKSO,'function_handle')
        fprintf('  * Output functions:\n'); 
        if ~isempty(obj.sFRho), fprintf('    > Output function of the density\n');          end
        if ~isempty(obj.sFKSO), fprintf('    > Output function of the orbitals\n');         end
    end

    if obj.sDFT,    fprintf('  * Save the DFT object in individual files\n');               end

    ref                 =   [ref,obj.docOutputChildren];

    fprintf('\n');
    
    % Bibliography, funding, and footer ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QMol_doc.showBibliography(ref);
    QMol_doc.showFunding;
    if ~obj.isRun,          QMol_doc.showFooter;                            end
    
end