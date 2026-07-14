%% make_workflow_schematics.m
%  Renders the two preprocessing workflow schematics as JPG images.
%  Pure plotting (no toolboxes, no data) -> run and collect the .jpg files.
%
%  Outputs (next to this script):
%    fiber_preprocessing_schematic.jpg
%    cellular_preprocessing_schematic.jpg
%
%  Design: each box's height is derived from its number of body lines (see
%  box()), so text can never overflow its box regardless of content changes.
%  Keep body text to short, one-line-per-step summaries -- this is a
%  workflow overview, not an implementation walkthrough (see each
%  pipeline's own README.md for full detail).

function make_workflow_schematics()
    out_dir = fileparts(mfilename('fullpath'));
    draw_fiber(fullfile(out_dir, 'fiber_preprocessing_schematic.jpg'));
    draw_cellular(fullfile(out_dir, 'cellular_preprocessing_schematic.jpg'));
    fprintf('Done. Wrote 2 JPGs to %s\n', out_dir);
end

% ---------- palette ----------
function c = col(name)
    switch name
        case 'input',   c = [0.86 0.92 0.98];  % light blue
        case 'imaging', c = [0.88 0.95 0.86];  % light green
        case 'ephys',   c = [0.98 0.91 0.80];  % light orange
        case 'proc',    c = [0.93 0.89 0.97];  % light purple
        case 'analysis',c = [0.99 0.93 0.93];  % light red
        case 'output',  c = [0.85 0.87 0.90];  % grey
        case 'edge',    c = [0.25 0.25 0.30];
        otherwise,      c = [0.95 0.95 0.95];
    end
end

% ---------- primitives ----------
function h = box(x, y_top, w, ttl, body_lines, fillname)
%BOX  Draw a box whose height is derived from its content, so text can
%     never overflow. body_lines is a cell array of short strings (one per
%     line); pass {} for a title-only box. Returns the box height so the
%     caller can stack the next box below it.
    TITLE_H  = 0.50;
    LINE_H   = 0.36;
    PAD_TOP  = 0.24;
    PAD_BOT  = 0.20;

    n = numel(body_lines);
    h = PAD_TOP + TITLE_H + n*LINE_H + PAD_BOT;
    y = y_top - h;

    rectangle('Position',[x y w h],'Curvature',0.12,'FaceColor',col(fillname), ...
        'EdgeColor',col('edge'),'LineWidth',1.3);
    text(x+0.20, y_top-PAD_TOP, ttl, 'FontWeight','bold','FontSize',10.5, ...
        'VerticalAlignment','top','Interpreter','none');
    for i = 1:n
        ty = y_top - PAD_TOP - TITLE_H - (i-1)*LINE_H;
        text(x+0.20, ty, body_lines{i}, 'FontSize',9, ...
            'VerticalAlignment','top','Interpreter','none','Color',[0.15 0.15 0.15]);
    end
end

function arrow(x, y1, y2)
    % vertical down arrow from y1 (top) to y2 (bottom)
    annline(x, y1, x, y2);
    hw = 0.10; hh = 0.16;
    patch([x-hw x+hw x], [y2+hh y2+hh y2], col('edge'), 'EdgeColor','none');
end

function annline(x1,y1,x2,y2)
    line([x1 x2],[y1 y2],'Color',col('edge'),'LineWidth',1.4);
end

function newfig(ttl, W, H)
    figure('Color','w','Units','pixels','Position',[60 60 W*95 H*95]);
    axes('Position',[0 0 1 1]); hold on; axis off;
    xlim([0 W]); ylim([0 H]);
    text(W/2, H-0.35, ttl, 'FontWeight','bold','FontSize',15, ...
        'HorizontalAlignment','center','Interpreter','none');
end

function save_jpg(fig_path, y_used, H)
    % Trim unused canvas below the last box before saving.
    ylim([max(y_used-0.3, 0), H]);
    set(gcf,'PaperPositionMode','auto');
    print(gcf, fig_path, '-djpeg', '-r150');
    close(gcf);
    fprintf('  wrote %s\n', fig_path);
end

% ============================================================== FIBER
function draw_fiber(fig_path)
    W = 12; H = 20;
    newfig('Fiber Photometry + Open Ephys -- Preprocessing Workflow', W, H);

    bw = 10.6; bx = 0.7;
    gap = 0.35;
    y = H-1.0;

    % --- inputs: two side-by-side ---
    iw = 5.15;
    h1 = box(bx, y, iw, 'Imaging', {'GEVI voltage, multi-site optic fibers'}, 'imaging');
    box(bx+iw+0.3, y, iw, 'Electrophysiology', {'LFP, stimulation, camera trigger, running speed'}, 'ephys');
    y = y - h1 - gap;

    arrow(bx+bw/2, y+gap, y);

    h = box(bx, y, bw, 'Configure & set up', {'Load run parameters; add external toolboxes to path'}, 'input');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Load imaging', {'Concatenate raw frames; optional motion correction'}, 'imaging');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Extract fiber traces', {'Manually defined per-fiber ROI -> mean fluorescence trace'}, 'imaging');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Remove artifacts', {'Clean electrical/stimulation artifacts from each trace'}, 'proc');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Align to electrophysiology', ...
        {'Detect stimulation onset; align LFP and fiber traces to a', ...
         'common imaging timebase via camera-trigger pulses'}, 'ephys');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Correct photobleaching & compute deltaF/F', ...
        {'Double-exponential decay fit; normalize to a baseline window'}, 'proc');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Spectral analysis', {'Band power and time-frequency spectra per fiber'}, 'analysis');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Fiber-LFP coupling', ...
        {'Phase-locking and coherence, split by rest/running state'}, 'analysis');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Generate QC figures', {'Trace, spectrogram, and coupling diagnostics'}, 'output');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Save output', ...
        {'One file per trial: aligned traces, LFP, and metadata'}, 'output');
    y = y-h;

    save_jpg(fig_path, y, H);
end

% ============================================================== CELLULAR
function draw_cellular(fig_path)
    W = 12; H = 17.5;
    newfig('Single-Cell Voltage Imaging (DBS) -- Preprocessing Workflow', W, H);

    bw = 10.6; bx = 0.7;
    gap = 0.35;
    y = H-1.0;

    iw = 5.15;
    h1 = box(bx, y, iw, 'Imaging (per trial)', {'Single-neuron GEVI, ROI drawn per cell'}, 'imaging');
    box(bx+iw+0.3, y, iw, 'Electrophysiology', {'LFP, stimulation, camera trigger, motion'}, 'ephys');
    y = y - h1 - gap;
    arrow(bx+bw/2, y+gap, y);

    h = box(bx, y, bw, 'Configure & detect trials', ...
        {'Load run parameters and session metadata; find trial folders'}, 'input');
    y = y-h; arrow(bx+bw/2, y, y-gap); y = y-gap;

    % --- per-trial loop frame: reserve a dedicated label band above the
    %     dashed border, rather than relying on clearance above the border
    %     line (which let the label text touch it) ---
    label_h = 0.55;
    text(bx-0.1, y-0.08, 'For each trial (ROIs shared from the first trial)', ...
        'FontWeight','bold','FontSize',10,'Color',[0.45 0.30 0.55], ...
        'VerticalAlignment','top','Interpreter','none');
    loop_top = y - label_h;
    yt = loop_top - 0.30;

    yt2 = box(bx+0.15, yt, bw-0.3, 'Load imaging & correct motion', ...
        {'Concatenate frames; ROIs drawn once and shared across trials'}, 'imaging');
    yt = yt - yt2; arrow(bx+bw/2, yt, yt-0.28); yt = yt-0.28;

    yt2 = box(bx+0.15, yt, bw-0.3, 'Extract per-neuron traces', ...
        {'Mean fluorescence within each neuron ROI'}, 'imaging');
    yt = yt - yt2; arrow(bx+bw/2, yt, yt-0.28); yt = yt-0.28;

    yt2 = box(bx+0.15, yt, bw-0.3, 'Align to electrophysiology', ...
        {'Detect stimulation onset; align LFP to the imaging timebase'}, 'ephys');
    yt = yt - yt2; arrow(bx+bw/2, yt, yt-0.28); yt = yt-0.28;

    yt2 = box(bx+0.15, yt, bw-0.3, 'Correct photobleaching', ...
        {'Double-exponential decay fit per neuron'}, 'proc');
    yt = yt - yt2; arrow(bx+bw/2, yt, yt-0.28); yt = yt-0.28;

    yt2 = box(bx+0.15, yt, bw-0.3, 'Detect spikes', ...
        {'Per-neuron spike times, firing rate, subthreshold signal'}, 'analysis');
    yt = yt - yt2; arrow(bx+bw/2, yt, yt-0.28); yt = yt-0.28;

    yt2 = box(bx+0.15, yt, bw-0.3, 'Store trial result', {}, 'output');
    yt = yt - yt2;

    loop_h = loop_top - yt + 0.25;
    rectangle('Position',[bx-0.25 (loop_top-loop_h) bw+0.5 loop_h], ...
        'Curvature',0.04,'EdgeColor',[0.45 0.30 0.55],'LineWidth',2.0,'LineStyle','--');

    y = loop_top - loop_h - 0.15;
    arrow(bx+bw/2, y, y-gap); y = y-gap;

    h = box(bx, y, bw, 'Save output', ...
        {'One file per session: shared ROIs + per-trial signals, spikes, metadata'}, 'output');
    y = y-h;

    save_jpg(fig_path, y, H);
end
