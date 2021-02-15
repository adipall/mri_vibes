function solve_second_order_system(p, dt, nt, plot_scale, cases, f_xt, do_movie)

plot_frequency = 1;

nc = length(cases);
for ic = 1:nc
    n(ic) = size(cases{ic}.K,1);
    xm{ic} = zeros(2*n(ic),1);
    xi{ic} = zeros(2*n(ic),1);

    if isfield(cases{ic},'v0_ind') && isfield(cases{ic},'v0_val')
        xm{ic}(n(ic) + cases{ic}.v0_ind) = cases{ic}.v0_val;
    end

    cases{ic}.Kx = [cases{ic}.K, sparse(n(ic),n(ic));
                    sparse(n(ic),n(ic)), -speye(n(ic))];
    cases{ic}.Mx = [sparse(n(ic),n(ic)), cases{ic}.M;
                    speye(n(ic)), sparse(n(ic),n(ic))];
    cases{ic}.Fx = [cases{ic}.V'*cases{ic}.F; zeros(n(ic), 1)];
    cases{ic}.Ax = eye(2*n(ic)) + dt*(cases{ic}.Mx\cases{ic}.Kx);
end

err = zeros(nc,1);

set(0, 'DefaultLineLineWidth', 6);
set(0,'defaultAxesFontSize',20)

if do_movie
    vidObj = VideoWriter('movie.avi');
    open(vidObj);
end

for ic=1:nc
    labels{ic} = cases{ic}.label;
end

figure(1), clf
for it = 1:nt
    ti = (it-1)*dt;

    fi = f_xt(p,ti);

    for ic=1:nc
        xi{ic} = cases{ic}.Ax\(xm{ic} + dt*(cases{ic}.Mx\cases{ic}.Fx)*...
            cases{ic}.force_time(ti));
        xm{ic} = xi{ic};
    end
    
    if mod(it,plot_frequency) == 0
        figure(1), clf
        for ic=1:nc
            u_plot = cases{ic}.disp_vec*cases{ic}.force_time(ti);
            u_plot(cases{ic}.free) = cases{ic}.V*xi{ic}(1:n(ic));

            err(ic) = max(err(ic), norm(u_plot - fi));

            plot(p, u_plot), hold on
        end

        plot(p, fi, 'k-.')

        ylim(plot_scale)
        legend(labels)
        drawnow
        pause(.05)
        if do_movie
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
        end
    end
end

if do_movie
    close(vidObj);
end

for ic=1:nc
    fprintf('Case %d (%s):\n',ic,labels{ic});
    fprintf('cond(M) = %g\n',cond(full(cases{ic}.M)));
    fprintf('max error: %g\n',err(ic))
end
