```
inputs: iterations, stepSize, stepSizeAngularMax, stepSizeAngularMin, eta, T, thetaMin, thetaMax, nBeams, beams[1..nBeams], W, d

struct
{
    float theta;
    float phi;
    float fluenceMap[][];
    float dose[][][];
} beam;

float doseLoss[];
float smoothnessLoss[];
float perturbationEnergy[];

// initial dose calculation
for b = 1,2,...nBeams do
    beams[b].fluenceMapToDose();
end for

for iter = 1,2,...iterations do
    // beam[b].dose is up-to-date for b = 1,2,...nBeams
    // calculate totalDose, doseLoss, and doseGrad
    totalDose := 0;
    for b = 1,2,...nBeams do
        totalDose += beams[b].dose;
    end for

    doseLoss[iter], doseGrad := calcDoseGrad(totalDose, W, d);

    for b = 1,2,...nBeams do
        fluenceGrad := beams[b].calcFluenceGrad(doseGrad);
        smoothnessLoss[(iter-1)*nBeams+b], smoothnessGrad := beams[b].calcSmoothnessGrad();
        totalGrad := fluenceGrad + eta * smoothnessGrad;
        totalGradNorm := norm(totalGrad);
        beams[b].fluenceMap -= step_size * totalGrad / totalGradNorm;

        // update beam.dose
        beams[b].fluenceMapToDose();
    end for

    // update T, as ln(2) times the average absolute difference of 
    // perturbation Energy of last iter
    if iter > 1 then
        absoluteDiff := 0;
        for b = 1,2,...nBeams do
            absoluteDiff += abs(perturbationEnergy[2*((iter-2)*nBeams+b)] - \
                perturbationEnergy[2*((iter-2)*nBeams+b)-1]);
        end for
        T = absoluteDiff / nBeams * ln(2);
    end if

    // update stepSizeAngular, as a linear decay
    stepSizeAngular = (stepSizeAngularMax * (iterations - iter) + \
        stepSizeAngularMin * iter) / iterations;
    
    // do perturbation
    for b = 1,2,...nBeams do
        // calculate loss for original beam angles
        totalDose := 0;
        for bb = 1,2,...nBeams do
            totalDose += beams[bb].dose;
        end for
        doseLoss0, _ := calcDoseGrad(totalDose, W, d);
        theta_org = beams[b].theta;
        phi_org = beams[b].phi;

        // update beam angles
        beams[b].theta += (random01() - 0.5) * stepSizeAngular;
        beams[b].phi += (random01() - 0.5) * stepSizeAngular;

        // clamp theta value
        if beams[b].theta < thetaMin then
            beams[b].theta = thetaMin;
        end if
        if beams[b].theta > thetaMax then
            beams[b].theta = thetaMax;
        end if
        beams[b].fluenceMapToDose();

        totalDose := 0;
        for bb = 1,2,...nBeams do
            totalDose += beams[bb].dose;
        end for
        doseLoss1, _ = calcDoseGrad(totalDose, W, d);

        if doseLoss1 <= doseLoss0 then
            probability = 1;
        else
            probability = exp((doseLoss0 - doseLoss1) / T);
        end if

        // if random01() <= probability, then take the new angles, do nothing
        if (random01() > probability) then
            // restore original angles
            beams[b].theta = theta_org;
            beams[b].phi = phi_org;
            beams[b].fluenceMapToDose();
        end if
    end for
end for
```