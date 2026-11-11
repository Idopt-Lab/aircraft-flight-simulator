function eng = SilverScytheEngineDeck(ss)
% SILVERSCYTHEENGINEDECK  NPSS engine deck for one F414-GE-400.
%
% Loads data/F414_NPSS_deck.csv, which was extracted from the team's NPSS
% run F414-GE-400_Outfile.txt (NPSS 3.2, model cdm01_tf02, 16 April 2026).
% The file holds 1,039 converged cases on a complete
%
%   13 Mach (0 to 0.60)  x  8 altitude (0 to 35,000 ft)  x  10 PLA
%
% grid, with net thrust, TSFC, corrected airflow and overall pressure
% ratio.  Report section 13.6 describes this as the higher-fidelity engine
% analysis used to verify the Mattingly deck, which had errors of up to
% 25%.  Its design point matches report Table 11 exactly (168.8 lbm/s,
% 13,700 lbf dry / 0.8552, 22,000 lbf wet / 1.8439).
%
% Throttle mapping used by SilverScytheThrustModel:
%   throttle 0.0 -> PLA  10  (idle)
%   throttle 0.5 -> PLA  99  (maximum dry / military)
%   throttle 1.0 -> PLA 150  (maximum afterburning)
%
% NPSS coverage stops at Mach 0.60 and 35,000 ft, so the deck cannot on its
% own reach the Mach 1.64 dash point.  SilverScytheThrustModel extends it
% with the Mattingly model scaled to match NPSS on the coverage boundary.

if nargin < 1 || isempty(ss), ss = SilverScytheData(); end

persistent cache cache_key
key = char(fullfile(ss.paths.data_dir,char(ss.engine.npss_deck_file)));
if ~isempty(cache) && strcmp(cache_key,key)
    eng = cache;
    return;
end

if ~isfile(key)
    error('SilverScytheEngineDeck:MissingDeck', ...
        'NPSS deck not found: %s',key);
end

T = readtable(key);
required = ["Mach","Altitude_ft","PLA","Fn_lbf","TSFC_lbm_per_lbf_hr"];
names = string(T.Properties.VariableNames);
for k = 1:numel(required)
    if ~any(strcmpi(names,required(k)))
        error('SilverScytheEngineDeck:BadDeck', ...
            'NPSS deck is missing column %s.',required(k));
    end
end

mach = unique(T.Mach);
alt  = unique(T.Altitude_ft);
pla  = unique(T.PLA);

n_expected = numel(mach)*numel(alt)*numel(pla);
Fn   = nan(numel(mach),numel(alt),numel(pla));
TSFC = nan(numel(mach),numel(alt),numel(pla));
Wair = nan(numel(mach),numel(alt),numel(pla));

for r = 1:height(T)
    i = find(mach == T.Mach(r),1);
    j = find(alt  == T.Altitude_ft(r),1);
    k = find(pla  == T.PLA(r),1);
    Fn(i,j,k)   = T.Fn_lbf(r);
    TSFC(i,j,k) = T.TSFC_lbm_per_lbf_hr(r);
    if any(strcmpi(names,"Wair_lbm_per_s"))
        Wair(i,j,k) = T.Wair_lbm_per_s(r);
    end
end

n_filled = sum(isfinite(Fn(:)));
if n_filled < n_expected
    % Fill any gap along PLA by monotone interpolation so that the
    % gridded interpolant stays well defined.
    for i = 1:numel(mach)
        for j = 1:numel(alt)
            col = squeeze(Fn(i,j,:));
            good = isfinite(col);
            if nnz(good) >= 2 && ~all(good)
                Fn(i,j,~good) = interp1(pla(good),col(good),pla(~good), ...
                    'linear','extrap');
                colT = squeeze(TSFC(i,j,:));
                TSFC(i,j,~good) = interp1(pla(good),colT(good), ...
                    pla(~good),'linear','extrap');
            end
        end
    end
end

eng = struct();
eng.source        = key;
eng.mach_grid     = mach;
eng.alt_ft_grid   = alt;
eng.pla_grid      = pla;
eng.Fn_lbf        = Fn;
eng.TSFC          = TSFC;
eng.Wair_lbm_per_s = Wair;
eng.n_cases       = height(T);
eng.mach_max      = max(mach);
eng.alt_max_ft    = max(alt);
eng.pla_idle      = min(pla);
eng.pla_mil       = 99;
eng.pla_max_ab    = max(pla);
eng.Fn_fn   = griddedInterpolant({mach,alt,pla},Fn,'linear','nearest');
eng.TSFC_fn = griddedInterpolant({mach,alt,pla},TSFC,'linear','nearest');
eng.method  = "NPSS 3.2 engine deck (report section 13.6, Table 11)";

cache = eng;
cache_key = key;

end
