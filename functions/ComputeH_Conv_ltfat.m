classdef ComputeH_Conv_ltfat
    % if only the elements of h is changed, you can use this code.
    % if L or M or a is changed, you must initialize.
    properties
        wft         % FFT of shifted window
        gft         % FFT of shifted dual window
        Prod_FwG    % Product of wft and G
        frame
    end

    methods
        % コンストラクタ
        function obj = ComputeH_Conv_ltfat(gt,gd, L, M)
            b=L/M;
            % obj.wft = fft(circshift(buffer(gt, L), -M/2));
            % obj.gft = fft(circshift(buffer(gd, L), -M/2));
            obj.wft = fft([gt(1:M/2+1); zeros(L-M,1); gt(end-M/2+2:end)]);
            obj.gft = fft([gd(1:M/2+1); zeros(L-M,1); gd(end-M/2+2:end)]);

            q = 0:M-1;
            GIdx = mod((1:L)' - b*q - 1, L) + 1; % L = length(h)
            G = obj.gft(GIdx);
            % G(:, 2:2:end) = -1 * G(:, 2:2:end); % 偶数列の符号反転

            obj.Prod_FwG = obj.wft .* G;

            obj.frame=[];
        end

        % ComputeH メソッド
        function H = ComputeH(obj, h, a, M)
            b=length(h)/M;
            % h に関する計算
            hft = fft(h);
            % Hadamard product and downsampling before IFFT
            H=zeros(M,length(h)/a,M);
            for m=0:M-1
                H(m+1,:,:) = squeeze(ifft(sum(reshape( ...
                    hft .* circshift(obj.Prod_FwG,[b*m m]), length(h)/a,a,M), 2)));
            end
            H=H/a;
        end


        function h = InvComputeH(obj, gd ,gt, H, a, M)
            % h に関する計算
            L=size(obj.Prod_FwG,1);
            b=L/M;
            if isempty(obj.frame)
                obj.frame=zeros(a,1);
                win_conv=cconv(gt.^2, gd.^2,L);
                for l=1:a
                    temp=sum(win_conv(l:a:end));
                    obj.frame(l,1)=temp;
                end
            end
            
            % h = zeros(a,L/a);
            h = zeros(L,1);
            for m = 0:M-1
                %%%%%%%%%%%%%%%
                % tmp = squeeze(sum(reshape(circshift(obj.Prod_FwG,[b*m m]),[],a,M),2))/a;
                % for j=1:a
                % h(j,:) = h(j,:) + reshape(sum(ifft(circshift(tmp,j-1,1) .* fft(squeeze(H(m+1,:,:))) ), 2),1,[]);
                % end
                % h = h + reshape(sum(ifft( squeeze(sum(reshape(circshift(obj.Prod_FwG,[b*m m]),[],a,M),2))/a .* fft(squeeze(H(m+1,:,:))) ), 2),[],1);
                %%%%%%%%%%%%%%%
        

                h = h + sum(ifft( circshift(obj.Prod_FwG,[b*m m]) .* fft(upsample(squeeze(H(m+1,:,:)),a)) ), 2);
            end
            h=real(h);
            h=reshape(reshape(h,a,[])./obj.frame,[],1)/M^2;
        end
    end
end
            