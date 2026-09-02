# Working agreement

## How to talk to me (가장 중요)

**항상 하이레벨로, 쉽게.** 기계장치가 아니라 **뜻**을 말할 것.

- 먼저 "그래서 뭐가 달라지는데"를 한 줄로. 그 다음에 근거.
- 숫자는 before → after 로. 표가 문장보다 낫다.
- 함수 이름, 파일 경로, 변수명 나열은 답이 아니다. 물어보면 그때 말한다.
- 코드 블록은 그게 **핵심 논증**일 때만. 보여주기용 금지.
- 모르면 모른다고. 확인 안 했으면 확인 안 했다고.
- 끝내면서 "이거 할까요?"로 되묻지 말 것. 이미 지시받은 건 그냥 한다.

나쁜 예: "`theta_of_mu`에 `compensation='polaron'` 인자를 추가하고 `_deep`,
`phase_boundaries`, `distribute`에 스레딩했습니다."

좋은 예: "산소 하나 빠지면 전자 두 개가 남아서 Ti³⁺ 두 개가 됩니다. 그 자리도
셈에 넣으니까 결함 농도가 0.954에서 0.066으로 떨어졌습니다. 전에는 Ti₂O₃보다
네 배 더 환원된 상태를 예측하고 있었던 겁니다."

## What this repo is

Ti–O 고체–기체 열역학. 논문에 실제로 쓰이는 **계산 두 개**와, 그걸 보여주는
웹페이지 하나.

### 1. Gas–solid equilibrium (`solidgas/activeset.py`)

주어진 기체 charge에서 평형 전환율, 안정한 Ti–O 상, 가장 가까운 환원상까지의
여유를 준다. 응축상은 active set 으로 들고 가므로 **없는 상은 정확히 0** 이고,
얼마나 멀리 있는지는 KKT reduced cost 로 따로 보고한다. 루타일이 살아남지
못하는 feed 는 닫힌형으로 나온다.

### 2. Surface-capacity bound (`solidgas/vacancy_inventory.py`)

적정은 숫자 **하나**만 준다 — 빠진 산소 총량. 어디서 빠졌는지는 안 알려준다.
루타일 격자상수와 노출 면적만으로 (110) bridging 줄이 담을 수 있는 양을
계산하고, 거기서 **부등식**을 낸다: 측정값 중 최대 몇 %가 표면일 수 있는가,
따라서 최소 몇 %가 아래에 있어야 하는가.

- 0.9 µm 입자: bridging 줄 전체 13.6, 재배열 시 6.8 µmol-O g⁻¹
- 측정 95 µmol-O g⁻¹ → 표면은 **최대 7.1 %**, 아래가 **최소 92.9 %**
- (1×2) 0.5 ML 결손이 **유일한** 문헌 구조 입력. DFT 에너지는 안 들어간다.

**표면 아래를 subsurface / bulk / extended defect 로 나누지 않는다.** 깊이를
분해하는 측정이 이 연구에 없다. 나눴던 코드(3-class stat mech, 몬테카를로,
반경 프로파일, 수송, CS 용해도)는 전부 제거했다 — 위 부등식을 그대로 재현하면서
측정이 확인할 수 없는 가정만 얹고 있었다.

## Rules

- 작업 브랜치: `claude/python-thermodynamics-engine-lk1v2g`. 다른 데 푸시 금지.
- 배포 방식: 브랜치 푸시 → main fast-forward. 게이트가 빨간 채로 main 을 올리지 않는다.
- 테스트는 0 skipped / 0 xfailed. skip 으로 넘기지 않는다.
- PR 은 명시적으로 요청받았을 때만 만든다.
- 레포에 들어가는 어떤 글에도 모델 이름을 쓰지 않는다 (커밋 메시지, 코드 주석, 문서).
- 파일 작업은 가능하면 bash (cat/sed/grep/heredoc) 로.
- **브라우저 엔진은 반드시 parity 게이트와 짝이다.** 페이지는 계산한다 —
  격자를 미리 구워두면 사용자가 실제로 가진 조성이나 BET 면적에 답할 수 없다.
  대신 `web/*.js` 엔진마다 파이썬 원본에 붙들어 매는 테스트가 있어야 한다.
  게이트 없는 두 번째 구현을 만들지 말 것.
- **논문이 인용하는 숫자는 페이지에서 다시 계산하지 않는다.** 커밋된 JSON 을
  읽는다. 사용자가 뭘 입력하든 논문 값과 어긋날 수 없게.
- 논문이 쓰지 않고 논문 데이터로 맞출 수도 없는 모델은 레포에 두지 않는다.

## Verification standard

숫자를 주장하려면 **독립 구현 세 개가 일치**해야 한다.

1. 파이썬 엔진 (`solidgas/`)
2. mpmath 80자리 오라클 (`scripts/oracle_tio.py` → `data/reference_*.json`)
3. 브라우저 미러 (`web/activeset.js`, `web/vacancy.js`)

오라클은 엔진을 호출하지 않는다. 안 그러면 검증이 아니라 복사다.

브라우저 미러는 물리가 아니라 **포팅**을 검증한다. 그래도 필요하다 — 페이지가
계산을 하는 이상 두 구현이 조용히 갈라지는 걸 막아야 한다. 예산: activeset 은
4 ulp (glibc exp 대 V8 exp), vacancy 는 **0 ulp** (지수함수가 없다).

표면 용량 쪽은 오라클이 없다 — 솔버가 없기 때문이다. 대신 게이트가 산술을
검증한다: 단위격자가 루타일 밀도 4.25 g/cm³ 와 19.22 Å² (110) 셀을 재현하는지,
부등식이 부등식인지, 그리고 **DFT 에너지가 흘러들어오지 않았는지**.

## Reproducing

```bash
python3 scripts/reproduce_paper.py    # paper_outputs/*.csv + site_data.json
python3 scripts/build_site.py         # docs/index.html
python3 -m pytest tests/ -q
PW_DIR=<dir with node_modules/playwright> python3 scripts/check_page.py
```

`paper_outputs/` 는 커밋되어 있고, 재생성 시 한 바이트라도 달라지면 게이트가
떨어진다.
