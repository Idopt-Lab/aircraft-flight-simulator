function deck = SilverScytheAeroDeck(ss,options)
% SILVERSCYTHEAERODECK  Build (or load) the high-fidelity aerodynamic deck.
%
% Runs the report's higher-fidelity aerodynamic chain once over a
% Mach/altitude/CL grid and stores it as fast interpolants:
%
%   CD_profile(Mach, altitude)  FRICTION   - SilverScytheFriction
%   e_span, k_induced           LIDRAG     - SilverScytheLidrag
%   CD_wave(Mach, CL)           Raymer     - SilverScytheWaveDrag
%
% Tabulating is what makes the chain usable inside a 6-DOF simulation: the
% FRICTION Newton solve and the wave-drag curve fits are far too slow to
% run per force evaluation, and the tabulated form is also continuous in
% Mach, which the raw piecewise curve fits are not.
%
% The deck is cached in results/SilverScytheAeroDeck.mat and reused when the
% stored grid and geometry signature match.
%
% Usage
%   deck = SilverScytheAeroDeck();                       % load or build
%   deck = SilverScytheAeroDeck(ss,struct('rebuild',true));

if nargin < 1 || isempty(ss), ss = SilverScytheData(); end
if nargin < 2 || isempty(options), options = struct(); end
if ~isfield(options,'rebuild')  || isempty(options.rebuild),  options.rebuild = false; end
if ~isfield(options,'verbose')  || isempty(options.verbose),  options.verbose = true;  end

cache_path = fullfile(ss.paths.results_dir,char(ss.deck.cache_file));
signature = deck_signature(ss);

if ~options.rebuild && isfile(cache_path)
    stored = load(cache_path,'deck');
    if isfield(stored,'deck') && isfield(stored.deck,'signature') && ...
            isequal(stored.deck.signature,signature)
        deck = rebuild_interpolants(stored.deck);
        if options.verbose
            fprintf('Silver Scythe aero deck loaded from cache (%s).\n', ...
                ss.deck.cache_file);
        end
        return;
    end
end

if options.verbose
    fprintf('\nBuilding Silver Scythe high-fidelity aero deck\n');
    fprintf('  FRICTION : %d Mach x %d altitude points\n', ...
        numel(ss.deck.mach_grid),numel(ss.deck.altitude_ft_grid));
    fprintf('  LIDRAG   : %d spanwise stations\n',ss.deck.spanload_stations);
    fprintf('  wave drag: %d Mach x %d CL points\n', ...
        numel(ss.deck.mach_grid),numel(ss.deck.CL_grid));
end

mach_grid = ss.deck.mach_grid(:);
alt_grid  = ss.deck.altitude_ft_grid(:);
CL_grid   = ss.deck.CL_grid(:);

%% Profile drag by FRICTION
fr = SilverScytheFriction(ss,mach_grid,alt_grid);
CD_profile = fr.CD_profile.';           % (mach, altitude)
CD_friction = fr.CD_friction.';
CD_form = fr.CD_form.';

%% Span efficiency by LIDRAG
span = SilverScytheSpanload(ss);
lid = SilverScytheLidrag(span.eta,span.ccl_ca,ss.AR);

%% Wave drag over Mach and CL
CD_wave = zeros(numel(mach_grid),numel(CL_grid));
for j = 1:numel(CL_grid)
    w = SilverScytheWaveDrag(ss,mach_grid,CL_grid(j));
    CD_wave(:,j) = w.CD_wave(:);
end
w_ref = SilverScytheWaveDrag(ss,mach_grid,0);

deck = struct();
deck.signature   = signature;
deck.built_on    = datetime('now');
deck.mach_grid   = mach_grid;
deck.alt_ft_grid = alt_grid;
deck.CL_grid     = CL_grid;
deck.CD_profile  = CD_profile;
deck.CD_friction = CD_friction;
deck.CD_form     = CD_form;
deck.CD_wave     = CD_wave;
deck.e_span      = lid.e_span;
deck.k_induced   = lid.k_induced;
deck.lidrag      = lid;
deck.spanload    = span;
deck.form_factor = fr.form_factor;
deck.component_names = ss.friction.component_names;
deck.M_DD_zero_lift  = w_ref.M_DD_zero_lift;
deck.quarter_chord_sweep_deg = w_ref.quarter_chord_sweep_deg;
deck.EWD = ss.wave.EWD;
deck.method = "FRICTION + LIDRAG + Raymer wave drag (report sections 11.4-11.5)";

if ~isfolder(ss.paths.results_dir)
    mkdir(ss.paths.results_dir);
end
save(cache_path,'deck');

deck = rebuild_interpolants(deck);

if options.verbose
    fprintf('  e_span   : %.5f  (k = %.6f)\n',deck.e_span,deck.k_induced);
    fprintf('  deck cached to %s\n\n',cache_path);
end

end

% -----------------------------------------------------------------------
function deck = rebuild_interpolants(deck)
% griddedInterpolant objects do not survive save/load reliably across
% releases, so they are always rebuilt from the stored grids.
deck.CD_profile_fn = griddedInterpolant({deck.mach_grid,deck.alt_ft_grid}, ...
    deck.CD_profile,'linear','nearest');
deck.CD_wave_fn = griddedInterpolant({deck.mach_grid,deck.CL_grid}, ...
    deck.CD_wave,'linear','nearest');
end

% -----------------------------------------------------------------------
function sig = deck_signature(ss)
sig = struct();
sig.mach_grid   = ss.deck.mach_grid(:).';
sig.alt_ft_grid = ss.deck.altitude_ft_grid(:).';
sig.CL_grid     = ss.deck.CL_grid(:).';
sig.stations    = ss.deck.spanload_stations;
sig.geometry    = [ss.S_ref_ft2 ss.span_ft ss.root_chord_ft ss.tip_chord_ft ...
                   ss.LE_sweep_deg ss.TE_sweep_deg ss.tc_wing ...
                   ss.max_cross_section_ft2 ss.wave_reference_length_ft];
sig.friction    = [ss.friction.S_wet_ft2(:).' ss.friction.l_ref_ft(:).' ...
                   ss.friction.t_c(:).' ss.friction.geom_type(:).' ...
                   ss.friction.Q(:).' ss.friction.transition(:).' ...
                   ss.friction.TwTaw ss.friction.leakage_multiplier];
sig.EWD         = ss.wave.EWD;
sig.version     = "v1";
end
