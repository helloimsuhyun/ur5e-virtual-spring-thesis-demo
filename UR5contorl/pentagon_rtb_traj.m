function x_ref = pentagon_rtb_traj(t)
% x_ref = pentagon_rtb_traj(t)
%   t      : 현재 시간 [s]
%   x_ref  : 3x1, [x; y; z]
%
% 내부에서 RTB의 mstraj를 한 번만 써서 오각형 궤적을 만들고,
% 이후에는 interp1로 시간 t에 맞는 x_ref(t)를 반환함.

    persistent tvec x_table T_total

    if isempty(tvec)
        % ----- 궤적 한 바퀴 설정 -----
        T_total = 10;   % 한 바퀴 도는 시간 [s]
        dt      = 0.01; % 궤적 샘플링 시간 [s]

        [tvec, x_table] = makePentagonTrajRTB(T_total, dt);
        % tvec : 1xN 시간 벡터
        % x_table : N×3, 각 행이 [x y z]
    end

    % 주기적으로 반복되도록 시간 랩핑
    tau = mod(t, T_total);

    % 선형 보간으로 현재 위치 계산
    xyz = interp1(tvec, x_table, tau, 'linear', 'extrap');  % 1x3

    x_ref = xyz(:);   % 3x1로 반환
end


% ========= 여기부터는 보조 함수 =========
function [tvec, x_table] = makePentagonTrajRTB(T_total, dt)
% [tvec, x_table] = makePentagonTrajRTB(T_total, dt)
%   T_total : 전체 한 바퀴 시간 (대충 맞춰줌, 실제로는 T_move+T_dwell로 결정)
%   dt      : 샘플링 시간 [s]
%   tvec    : 1×N 시간 벡터
%   x_table : N×3, 각 행이 [x y z]

    % ----- 오각형 파라미터 -----
    N = 5;                        % 꼭짓점 개수
    center = [0.5 0.0 0.4];       % 오각형 중심 [x y z]
    R = 0.30;                     % 크기(반지름)

    % ----- "한 변 이동 시간" / "꼭짓점에서 멈추는 시간" 설정 -----
    T_move  = 1.0;   % 각 변을 따라 이동하는 시간 [s] → 작을수록 빠르게
    T_dwell = 4.0;   % 꼭짓점에서 머무는 시간 [s]

    % 총 시간은 대략:
    T_total = N*(T_move + T_dwell);

    % ----- 기본 꼭짓점 (한 번씩) 계산: x 고정, y-z 평면에서 오각형 -----
    V = zeros(N,3);  % 꼭짓점들

    for k = 1:N
        theta = deg2rad(90 + (k-1)*360/N);   % 90도에서 시작, CCW

        px = center(1);                      % ✅ x 고정
        py = center(2) + R*cos(theta);       % y 오각형
        pz = center(3) + R*sin(theta);       % z 오각형

        V(k,:) = [px py pz];
    end

    % 마지막 꼭짓점에서 다시 첫 번째로 돌아오도록 (폐곡선)
    V = [V; V(1,:)];   % (N+1)×3

    % ----- via-point 확장: [이동 → 멈춤] 패턴 만들기 -----
    %   pvia:  시작점, (이동 끝점, 멈춤점)을 반복
    %   예) V1 → V2(이동), V2(멈춤), V3(이동), V3(멈춤), ...

    pvia = [];
    tseg = [];  % 각 구간 시간들

    % 시작점
    pvia = [pvia; V(1,:)];

    for k = 1:N
        % 이동 구간: V(k) → V(k+1)
        pvia = [pvia; V(k+1,:)];  % 이동 끝점 추가
        tseg = [tseg, T_move];

        % 멈춤 구간: V(k+1) → V(k+1) (같은 점 두 번)
        pvia = [pvia; V(k+1,:)];  % 같은 점 한 번 더
        tseg = [tseg, T_dwell];
    end
    % pvia 크기: (2N+1)×3
    % tseg 크기: 1×(2N)  (구간 개수 = 점 개수 - 1)

    % ----- RTB mstraj로 다구간 궤적 생성 -----
    % qdmax는 []로 두고, tseg로 각 구간 시간을 직접 지정
    qdmax   = [];                     % 속도 자동 (tseg에 의해 결정)
    q0      = pvia(1,:);              % 시작점
    tacc    = 0.3 * T_move;           % 가속 시간(적당히)

    x_table = mstraj(pvia, qdmax, tseg, q0, dt, tacc);

    % 시간 벡터 생성
    tvec = (0:size(x_table,1)-1) * dt;

    % 필요하면 T_total에 맞춰 잘라주기 (여기선 거의 안 써도 됨)
    if tvec(end) > T_total
        idx = tvec <= T_total;
        tvec    = tvec(idx);
        x_table = x_table(idx,:);
    end
end

