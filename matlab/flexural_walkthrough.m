%% Setup
clear; clc; close all

root   = 'Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons\P01-LT150-LH4.5';
cid    = 'P01-FCL00-02';
tag    = 'FCL0002';

dic_channel_qc(tag); % generate the DIC channel quality control figs

mtsdir = fullfile(root,'MTS');
rawdir = fullfile(root,'DIC','raw','2026_FSR_Flexural_FCL_FIS',tag);
tmpdir = fullfile(root,'DIC','tmp_flexural');
matdir = fullfile(tmpdir,'matlab');

geom   = readtable(fullfile(tmpdir,[cid '_geometry.csv']));
frames = readtable(fullfile(tmpdir,[cid '_frames.csv']));
levels = readtable(fullfile(matdir,[cid '_levels.csv']));

L = 8.00*25.4;
b = geom.b_mm;
d = geom.d_mm;

nlev = height(levels);
fld  = cell(nlev,1);
for k = 1:nlev
    fld{k} = readtable(fullfile(matdir,sprintf('%s_field_%d.csv',cid,k)));
end
lab = arrayfun(@(f) sprintf('%.2f kN',f/1000), levels.force_N, 'UniformOutput',false);

cmap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1)
        ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
cols = parula(nlev+1);

%% Section 1 - raw channels and their noise
mts   = readmatrix(fullfile(mtsdir,[cid '.txt']),'FileType','text','Delimiter','\t','NumHeaderLines',8);
i_mts = find(mts(:,2) == min(mts(:,2)),1);
t_mts = mts(1:i_mts,4);
f_mts = mts(1:i_mts,2);
d_mts = mts(1:i_mts,1);

sync  = readtable(fullfile(rawdir,[tag '.csv']),'VariableNamingRule','preserve');
i_dic = geom.break_frame + 1;
t_dic = sync.("Time_0_0")(1:i_dic) - sync.("Time_0_0")(1);
f_dic = sync.("Dev1/ai2")(1:i_dic);
d_dic = sync.("Dev1/ai1")(1:i_dic);

chan = {'MTS force',            'N',  t_mts, f_mts
        'MTS displacement',     'mm', t_mts, d_mts
        'DIC ai2 (nom. load)',  'V',  t_dic, f_dic
        'DIC ai1 (nom. disp)',  'V',  t_dic, d_dic};

fprintf('\n%-20s %6s %8s %10s %10s %10s %8s %8s\n', ...
        'channel','units','rate Hz','span','sigma_hf','sigma_res','ratio','SNR dB');
fprintf('%s\n',repmat('-',1,86));

figure('Name','Raw channels','Position',[80 80 1100 800])
for k = 1:4
    t = chan{k,3};
    x = chan{k,4};

    dt    = median(diff(t));
    win   = max(5,round(4/dt));
    trend = smoothdata(x,'sgolay',win);
    res   = x - trend;

    span  = max(trend) - min(trend);
    s_hf  = std(x(1:end-2) - 2*x(2:end-1) + x(3:end))/sqrt(6);
    s_res = std(res);

    fprintf('%-20s %6s %8.1f %10.4g %10.4g %10.4g %8.2f %8.1f\n', ...
            chan{k,1}, chan{k,2}, 1/dt, span, s_hf, s_res, s_res/s_hf, ...
            20*log10(span/s_hf));

    subplot(4,2,2*k-1)
    plot(t,x,'LineWidth',0.8)
    ylabel([chan{k,1} ' (' chan{k,2} ')'])
    grid on
    if k == 1, title('raw'), end
    if k == 4, xlabel('time (s)'), end

    subplot(4,2,2*k)
    plot(t,res,'LineWidth',0.6)
    ylabel(['residual (' chan{k,2} ')'])
    grid on
    if k == 1, title('detrended'), end
    if k == 4, xlabel('time (s)'), end
end

%% Section 2 - axial strain field at five load levels
lim = 0;
for k = 1:nlev
    v   = sort(abs(fld{k}.exx));
    lim = max(lim, v(round(0.99*numel(v)))*100);
end

figure('Name','Axial strain field','Position',[80 80 1000 850])
tiledlayout(nlev,1,'TileSpacing','compact','Padding','compact')
for k = 1:nlev
    nexttile
    scatter(fld{k}.X, fld{k}.Y, 6, fld{k}.exx*100, 'filled')
    colormap(cmap)
    clim([-lim lim])
    xline(geom.x_left_mm,'k:')
    xline(geom.x_right_mm,'k:')
    ylabel('Y (mm)')
    title(lab{k})
    if k == nlev, xlabel('position along specimen, X (mm)'), end
end
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = '\epsilon_{xx} (%)';

%% Section 3 - deflected shape
figure('Name','Deflected shape','Position',[80 80 800 550])
hold on
for k = 1:nlev
    T = fld{k};
    [g,xb] = findgroups(round(T.X/4)*4);
    n = splitapply(@numel,T.V,g);
    v = splitapply(@mean,T.V,g);
    xb = xb(n >= 5);
    v  = v(n >= 5);

    vl = interp1(xb,v,geom.x_left_mm,'linear','extrap');
    vr = interp1(xb,v,geom.x_right_mm,'linear','extrap');
    chord = vl + (xb - geom.x_left_mm)./(geom.x_right_mm - geom.x_left_mm).*(vr - vl);

    plot(xb, geom.defl_sign*(chord - v), 'LineWidth',1.6, 'Color',cols(k,:))
end
xline(geom.x_left_mm,'r:')
xline(geom.x_right_mm,'r:')
xlabel('position from center of span, X (mm)')
ylabel('deflection (mm)')
legend(lab,'Location','northwest')
grid on
hold off

%% Section 4 - curvature along the span
figure('Name','Curvature','Position',[80 80 800 550])
hold on
for k = 1:nlev
    T = fld{k};
    T = T(T.Y >= geom.y_roi_lo_mm + 0.30 & T.Y <= geom.y_roi_hi_mm - 0.30, :);

    [g,xb] = findgroups(round(T.X/5)*5);
    kap = nan(size(xb));
    for j = 1:numel(xb)
        m = (g == j);
        if sum(m) < 20, continue, end
        p = polyfit(T.Y(m), T.exx(m), 1);
        kap(j) = -p(1);
    end

    plot(xb, kap*1000, 'o', 'MarkerSize',4, 'MarkerFaceColor',cols(k,:), ...
         'MarkerEdgeColor','none')
end
yline(0,'k-')
xline(geom.x_left_mm,'r:')
xline(geom.x_right_mm,'r:')
xlabel('position from center of span, X (mm)')
ylabel('curvature \kappa (10^{-3} mm^{-1})')
legend(lab,'Location','northwest')
grid on
hold off

%% Section 5 - flexural stress-strain
keep = isfinite(frames.stress_MPa);
sig  = 3*frames.force_N(keep)*L./(2*b*d^2);

figure('Name','Stress-strain','Position',[80 80 800 600])
hold on
plot(frames.eps_curvature(keep)*100,  sig, 'LineWidth',1.6)
plot(frames.eps_deflection(keep)*100, sig, 'LineWidth',1.6)
plot(frames.eps_crosshead(keep)*100,  sig, 'LineWidth',1.6)
xlabel('flexural strain (%)')
ylabel('flexural stress (MPa)')
legend('DIC curvature','DIC deflection','MTS crosshead','Location','southeast')
xlim([0 inf])
ylim([0 inf])
grid on
hold off

fprintf('\n%s   sigma_fM = %.1f MPa   at %.2f %% strain\n', ...
        cid, max(sig), max(frames.eps_curvature(keep))*100);

%% Section 6 - curvature fit window sensitivity
% FlexuralDIC_Level1 reports one kappa per frame from a SINGLE exx-vs-Y fit
% pooling every point within MIDSPAN_HALF_WIDTH_MM = 10 mm of midspan. But
% kappa(x) is a triangle (Section 4), so a pooled fit returns the MEAN
% curvature over the window, not the apex value at midspan. For a triangle
% with zero crossings at +/-a,
%
%     kappa(w) / kappa(0) = 1 - w/(2a)
%
% so at w = 10 mm against a ~ 100 mm the reported kappa should sit ~5 % low,
% which inflates Ef_curvature by the same ~5 % and eats into the 8-11 % gap
% against Eq. (5). This sweeps w to measure the effect instead of assuming it.
% If the points fall on the dashed line, the artifact is confirmed and the
% correction is just the reciprocal of the formula above.

a_tri = 100;                        % kappa(x) zero crossing, mm (from Section 4)
xmid  = (geom.x_left_mm + geom.x_right_mm)/2;
hw    = [2.5 4 5 7.5 10 15 20 25];
[~,i10] = min(abs(hw - 10));
kap_hw  = nan(nlev,numel(hw));

for k = 1:nlev
    T = fld{k};
    T = T(T.Y >= geom.y_roi_lo_mm + 0.30 & T.Y <= geom.y_roi_hi_mm - 0.30, :);
    for j = 1:numel(hw)
        m = abs(T.X - xmid) <= hw(j);
        if sum(m) < 20, continue, end
        p = polyfit(T.Y(m), T.exx(m), 1);
        kap_hw(k,j) = -p(1);
    end
end

figure('Name','Curvature fit window','Position',[80 80 800 550])
hold on
for k = 1:nlev
    plot(hw, kap_hw(k,:)/kap_hw(k,1), 'o-', 'LineWidth',1.4, ...
         'Color',cols(k,:), 'MarkerFaceColor',cols(k,:))
end
plot(hw, (1 - hw/(2*a_tri))/(1 - hw(1)/(2*a_tri)), 'k--', 'LineWidth',1.2)
xline(hw(i10),'r:')
xlabel('fit half-width about midspan, w (mm)')
ylabel(sprintf('\kappa(w) / \kappa(%.1f mm)',hw(1)))
legend([lab; {sprintf('triangle, a = %g mm',a_tri)}],'Location','southwest')
grid on
hold off

fprintf('\n%-10s %12s %12s %10s %12s\n', ...
        'level','k(w_min)','k(10 mm)','ratio','Ef_c bias');
for k = 1:nlev
    r = kap_hw(k,i10)/kap_hw(k,1);
    fprintf('%-10s %12.4g %12.4g %10.4f %11.1f %%\n', ...
            lab{k}, kap_hw(k,1), kap_hw(k,i10), r, 100*(1/r - 1));
end
