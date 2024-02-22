function [test_mse,eva_mse]=ParticleFilter(N,steps)
    
    %initial
    last_velocity=1;
    last_heading=pi/4;
    heading_var=1/100;
    process_var=1;
    measurement_var=[1,1];
    xs=zeros(steps,3);
    xs_var=zeros(steps,3);
    x_initial=[0,0];
    %create particles
    particles =create_gaussian_particles(x_initial,last_heading, process_var, N);
    weights=ones([N,1])/N;
    %observe the dog move
    [z,position]=DogSimulation2Dcomplex(x_initial, last_velocity,heading_var,process_var,measurement_var,steps);
    %filter
    for i=1:steps

        %predict
        particles=predict(particles,heading_var,process_var,last_heading,last_velocity);
        %update
        weights=update(particles,weights,z(i,1:2),measurement_var);
        %resample
        if(neff(weights)<N/2)
            [particles,weights]=resample_systematic(particles,weights);
        end
        %estimate
        [mu,var]=estimate(particles,weights);
        last_heading=mu(3);
        if(i>1)
            diff_loc=mu(1:2)-xs(i-1,1:2);
        else
            diff_loc=mu(1:2);
        end
        last_velocity=sqrt(sum(diff_loc.^2));
        xs(i,:)=mu;
        xs_var(i,:)=var;
    end
    test_error = xs-position;
    test_mse =  sum(test_error.^2) / length(test_error );
    eva_error = z-position;
    eva_mse =  sum(eva_error.^2) / length(eva_error );
    plot(xs(:,1),xs(:,2));
    hold on
    plot(position(:,1),position(:,2));
    plot(z(:,1),z(:,2));
    legend("filter","track","measurement");
    hold off
    % plot(xs(:,3));
    % hold on
    % %plot(position(:,3));
    % plot(position(:,3));
    % legend("filter","track");
    % hold off
end
function particles=create_gaussian_particles(x,heading, initial_var, N)
    particles = [zeros([N, 2]),ones([N, 1])*heading];
    std=sqrt(initial_var);
    particles(:, 1) =x(1)+randn([N, 1]) * std;
    particles(:, 2) =x(2)+randn([N, 1]) * std;
end
function [mu,var]=estimate(particles,weights)
    %input:
    %particles: N*k
    %weight: N*1
    %output:
    %mu:1*k
    %var:1*k
    mu=sum(particles.*weights)/sum(weights);
    var=(particles-mu).^2;
    var=sum(var.*weights)/sum(weights);
end
function [particles]=predict(particles,heading_var,process_var,last_heading,last_velocity)
    %move according to control input u (heading change, velocity)   from
    %last move with noise Q (std heading change, std velocity)
    %particles (1:N,0): x position of N particles
    %particles (1:N,1): y position of N particles
    %particles (1:N,2): heading angle of N particles
    N=length(particles);
    %update heading
    particles(:, 3)=mod(last_heading+randn([N,1])*sqrt(heading_var),2*pi);
    % move in the (noisy) commanded direction |dt=1
    dist = last_velocity + randn([N,1]) * sqrt(process_var);
    % dist(dist>2)=2;
    % dist(dist<0)=0;
    particles(:, 1) =particles(:, 1)+cos(particles(:, 3)) .* dist;
    particles(:, 2) =particles(:, 2)+sin(particles(:, 3)) .* dist;
end
function [weights]=update(particles,weights,z,R)
    %get distances from particles
    %weights =weights.*(normpdf(z,particles(:,1),R(1))+normpdf(z,particles(:,2),R(2)));
    weights =weights.*(mvnpdf(z,particles(:,1:2),R));%2D gaussian
    weights =weights/ sum(weights);
end
function [particles,weights]=resample_systematic(particles,weights)
    N = length(weights);
    %make N subdivisions, choose positions with a random offset
    positions = (unifrnd(-1,0,[N,1])+(1:N)')/N;
    cumulative_sum = cumsum(weights);
    j=1;
    for i=1:N
        while(positions(i)>cumulative_sum(j))
            j=j+1;
        end
        particles(i,:)=particles(j,:);
    end
    %weights are equal after resampling
    weights=ones([N,1])/N;
end
function val=neff(weights)
    %when val<N/2, particles need resample
    val=1. /sum(weights.^2);
end