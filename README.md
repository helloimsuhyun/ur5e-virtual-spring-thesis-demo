# UR5e Virtual Spring-Damper Demo

본 프로젝트는 **UR5e 매니퓰레이터에 가상 스프링-댐퍼 모델을 적용한 제어 데모**이다.  
학부 전공 과목 **Robot Manipulation** 수업에서 학습한  
Virtual Spring-Damper control 방법을 MATLAB 및 Simscape 환경에서 구현하였다.

이 데모는 이론에서 다루는 dynamics를
시뮬레이션 상에서 직관적으로 확인하고 , 이해를 돕기 위한 예제이다.

---

## 📷 데모 화면
<p align="center">
  <img src="./example_gif.gif" width="745"/>
</p>

---

## 주요 기능

- Cartesian space에서의 virtual spring–damper 기반 position control
- Orientation error를 이용한 orientation control
- Simscape Multibody 기반 UR5e motion visualization
- 단순한 MATLAB code 구조로 수정 및 확장 용이

---

## 실행 방법

```bash
git clone https://github.com/<your-id>/ur5e-virtual-spring-damper-demo.git
cd ur5e-virtual-spring-damper-demo/UR5control
