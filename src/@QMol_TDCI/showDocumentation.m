function showDocumentation(obj)
%showDocumentation displays the run-time documentation
%   The run-time documentation prints out the specific configuration of the
%   TDCI object, as defined by its editable properties. The TDCI object
%   must be initialized before using showDocumentation.
%
%   obj.showDocumentation;
    
    % Header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ~obj.isRun
        QMol_doc.showHeader;
    end

    % Time propagation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QMol_doc.showSection('Time propagation');
    fprintf('  * The TDCI model dynamics relies on a canonical Hamiltonian formalism\n');
    fprintf('    [Mauger 2024, Mauger 2026].\n')

    ref                 =   [{'Mauger 2024','Mauger 2026'},obj.showDoc];

    fprintf('  * Parameters\n');
    fprintf('    Time interval = %s to %s a.u. (%s to %s fs)\n',...
        num2str(obj.T(1),'%5.2f'),                              ...
        num2str(obj.T(end),'%5.2f'),                            ...
        num2str(convertUnit.au2second(obj.T(1))*1e15,'%5.2f'),  ...
        num2str(convertUnit.au2second(obj.T(end))*1e15,'%5.2f'));
    fprintf('    Time step     = %s a.u. (%s as)\n',          ...
        num2str(obj.dt,'%5.3f'),                                ...
        num2str(convertUnit.au2second(obj.dt)*1e18,'%5.3f'));

    if ~isempty(obj.EF),    ref =   [ref,obj.EF.showDocumentation];                     end

    if ~isempty(obj.ABC),   ref =   [ref,obj.ABC.showDocumentation];                    end

    if obj.sRest,       fprintf('  * Restart option is enabled\n');
    else,               fprintf('  * Restart option is disabled\n');                    end

    fprintf('\n');
    
    % Output results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QMol_doc.showSection('Output results');

    if obj.sECI   ||   obj.sEWfcn
        fprintf('  * Configuration interaction and wave-function energies:\n'); 
        if obj.sECI,    fprintf('    > Save the configuration interaction energy\n');   end
        if obj.sEWfcn,  fprintf('    > Save the wave functions'' energy\n');            end
    end

    if obj.sDip   ||   obj.sVel   ||   obj.sAcc
        fprintf('  * Dipole, dipole velocity, and dipole acceleration:\n'); 
        if obj.sDip,    fprintf('    > Save the dipole signal\n');                      end
        if obj.sVel,    fprintf('    > Save the dipole-velocity signal\n');             end
        if obj.sAcc,    fprintf('    > Save the dipole-acceleration signal\n');         end
    end

    if obj.sIon
        fprintf('  * Ionization statistics:\n    > Save the ionization signal\n');  
    end

    if obj.sWfcn
        fprintf('  * Wave functions:\n'); 
        if obj.sWfcn,   fprintf('    > Save the wave functions\n');                     end
    end

    if isa(obj.sF,'function_handle')
        fprintf('  * Output function\n    > Output function of the wave functions\n');
    end

    if obj.sCI,     fprintf('  * Save the CI object in individual files\n');            end

    ref                 =   [ref,obj.docOutputChildren];

    fprintf('\n');
    
    % Bibliography, funding, and footer ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QMol_doc.showBibliography(ref);
    QMol_doc.showFunding;
    if ~obj.isRun,          QMol_doc.showFooter;                            end
    
end