function out = CalcWVDistances(cfg_in, out, out2)
% out = CalcWVDistances(cfg_in, out, out2)
%
% compute distances between all cell waveforms in out and all cell waveforms in out2
% (out and out2 are the output of CategorizeStriatumWave.m)
%
cfg_def = [];
cfg_def.method = 'SSE'; % {'SSE', 'SSEn', 'abs', 'absn', 'corr'}

cfg = ProcessConfig(cfg_def, cfg_in);


n1 = size(out.wv, 1);
n2 = size(out2.wv, 1);


for iX = 1:n1
    for iY = 1:n2
        
        wv1 = out.wv(iX, :);
        wv2 = out2.wv(iY, :);
        
        switch cfg.method
            case 'SSE'
                out.wv_dist(iX, iY) = sum((wv1 - wv2).^2);
            case 'SSEn' % normalize by energy of waveform 1
                out.wv_dist(iX, iY) = sum((wv1 - wv2).^2);
                out.wv_dist(iX, iY) = out.wv_dist(iX, iY) ./ sum((wv1 + wv2).^2);
            case 'abs'
                out.wv_dist(iX, iY) = sum(abs(wv1 - wv2));
            case 'absn'
                out.wv_dist(iX, iY) = sum(abs(wv1 - wv2));
                out.wv_dist(iX, iY) = out.wv_dist(iX, iY) ./ sum((wv1 + wv2));
            case 'corr'
                temp = corrcoef(wv1, wv2);
                out.wv_dist(iX, iY) = temp(1, 2);
        end
        
    end
end
