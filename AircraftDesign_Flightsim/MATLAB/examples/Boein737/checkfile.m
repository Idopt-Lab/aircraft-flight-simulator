requested_mach = [0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.50 0.60 0.65 0.70 0.72 0.76 0.82];

expected_usable_mach = [0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.50 0.60];

loaded_mach = data.grid.mach_vec(:).';

assert(numel(data.source_files) == 14, 'The parser did not receive all 14 files.');

assert(numel(loaded_mach) == numel(expected_usable_mach) && all(abs(loaded_mach-expected_usable_mach) < 1e-8), 'Usable DATCOM Mach grid is unexpected.');

missing_mach = requested_mach(~arrayfun(@(m) any(abs(loaded_mach-m) < 1e-8),requested_mach));

fprintf('\n=== B737 DATCOM CHECK ===\n');
fprintf('Files read          : %d\n',numel(data.source_files));
fprintf('Usable Mach values  :');
fprintf(' %.2f',loaded_mach);
fprintf('\n');
fprintf('Unavailable stations:');
fprintf(' %.2f',missing_mach);
fprintf('\n');
fprintf('Static tables       : %d\n',numel(data.static_tables));
fprintf('Dynamic tables      : %d\n',numel(data.dynamic_tables));

assert(numel(data.static_tables) == 90, 'Expected 90 usable static tables.');

assert(numel(data.dynamic_tables) == 90, 'Expected 90 usable dynamic tables.');

fprintf('Current valid DATCOM database check passed.\n');