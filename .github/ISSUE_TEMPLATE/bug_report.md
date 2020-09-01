---
name: Bug report
about: Create a report to fix bugs
title: ''
labels: bug
assignees: ''

---

## Description
Provide a clear and concise description of what the bug is.
Also give a description of how to fix the bug.


### To Reproduce:
What compilers/machines are you seeing this with?
Give explicit steps to reproduce the behavior.
1. do this
2. then that
3. then, oops, look at the bug


## Additional context
Add any other context about the problem here.
Directly reference any issues or PRs in this or other repositories that this is related to, and describe how they are related. Example:
- needs to be fixed also in noaa-emc/nems/issues/<issue_number>
- needed for noaa-emc/fv3atm/pull/<pr_number>


## Output

**Screenshots**
If applicable, drag and drop screenshots to help explain your problem.

**output logs**
If applicable, include relevant output logs.
Either drag and drop the entire log file here (if a long log) or

```
paste the code here (if a short section of log)
```

## Testing:

1. Have you tested the code changes? On what platforms?

2. Have you run regression test in ufs-weather-model or ufs-s2s-model with code changes?
- Will the baseline results change? 
- If the baseline results change, is it expected? Please give brief explanation.

## Dependent PRs:

Directly reference any issues or PRs in this or other repositories that this is related to, and describe how they are related. Example:
- required to support noaa-emc/GFDL_atmos_cubed_sphere/issues/<issue_number>
- ncar/ccpp-physics/pull/<pr_number>
- associated ufs-weather-model/pull/<pr_number>
