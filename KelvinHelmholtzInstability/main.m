//
//  main.m
//  KelvinHelmholtzInstability
//
//  Created by Jeffrey J. Early on 8/2/12.
//
//	This code is inspired by chapter 10.4 and problem 16.5 in Cushman-Roisin & Beckers "Introduction to Geophysical Fluid Dynamics", second edition.
//	I use the same parameters that they use in the 'shearedflow.m' example that can be found online.

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

int main (int argc, const char * argv[])
{
    @autoreleasepool {
        /************************************************************************************************/
        /*		Adjustable parameters																	*/
        /************************************************************************************************/
        
        GLFloat halfHeight = 10; // half the domain height, in meters
        GLFloat halfShearHeight = 0.1*halfHeight; // half the shear layer width, in meters
        GLFloat halfWidth = 4*halfHeight;
        GLFloat U0 = 1; // maximum fluid speed, in meters per second
        GLFloat dampingCoefficient = 0.05; // damping coefficient in meters per second, gets multiplied by the grid size.
        NSUInteger nPoints = 128; // Number of points in the height
        GLFloat maxTime = 200; // seconds
        GLFloat outputInterval = 1; // seconds
        
        /************************************************************************************************/
        /*		Define the problem dimensions															*/
        /************************************************************************************************/
        
        GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 4*nPoints domainMin: -halfWidth length: halfWidth*2.0];
        xDim.name = @"x";
        xDim.units = @"meters";
        GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLInteriorGrid nPoints: nPoints domainMin: -halfHeight length: halfHeight*2.0];
        yDim.name = @"y";
        yDim.units = @"meters";
        GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 1 + round(maxTime/outputInterval)  domainMin:0 length:maxTime];
        tDim.name = @"t";
        tDim.units = @"seconds";
        
        // Variables are always tied to a particular equation---so we create an equation object first.
        GLEquation *equation = [[GLEquation alloc] init];
        
        // Create variables associated with the two spatial dimensions
        NSArray *spatialDimensions = @[xDim, yDim];
        GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
        GLFunction *y = [GLFunction functionOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
        
        /************************************************************************************************/
        /*		Create and cache the differential operators that we will be using.						*/
        /************************************************************************************************/
        
        NSArray *differentiationBasis = @[@(kGLExponentialBasis),@(kGLSineBasis)];
        NSArray *spectralDimensions = [x dimensionsTransformedToBasis: differentiationBasis];
                
        // Create the operator xx+yy---this is how you compute y from eta
        GLLinearTransform *laplacian = [GLLinearTransform harmonicOperatorFromDimensions: spectralDimensions forEquation: equation];
        
        // Create the operator 1/(xx+yy)---this is how you compute eta from y.
        // Because we're in a *sine* basis, there is no zero wavenumber and this inversion works without dividing by zero.
        GLLinearTransform *inverseLaplacian = [laplacian inverse];
        
        GLLinearTransform *diff_xxx = [GLLinearTransform differentialOperatorWithDerivatives:@[@(3),@(0)] fromDimensions:spectralDimensions forEquation:equation];
        GLLinearTransform *diff_xyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(1),@(2)] fromDimensions:spectralDimensions forEquation:equation];
        GLLinearTransform *diff_xxy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(2),@(1)] fromDimensions:spectralDimensions forEquation:equation];
        GLLinearTransform *diff_yyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(0),@(3)] fromDimensions:spectralDimensions forEquation:equation];
        
        GLLinearTransform *diffJacobianX = [diff_xxx plus: diff_xyy];
        GLLinearTransform *diffJacobianY = [diff_xxy plus: diff_yyy];
        
        GLFloat k = 0.5*xDim.sampleInterval;
        GLLinearTransform *biharmonic = [GLLinearTransform harmonicOperatorOfOrder: 2 fromDimensions: spectralDimensions forEquation: equation];
        GLLinearTransform *svv = [GLLinearTransform spectralVanishingViscosityFilterWithDimensions: spectralDimensions scaledForAntialiasing: YES forEquation: equation];
        GLLinearTransform *diffLin = [[biharmonic times: @(k)] times: svv];
        
        /************************************************************************************************/
        /*		Create the initial shear flow															*/
        /************************************************************************************************/
        
        GLFunction *initialShearFlow = [GLFunction functionOfRealTypeWithDimensions: @[xDim, yDim] forEquation: equation];
        [initialShearFlow zero];
        for (NSUInteger i=0; i<initialShearFlow.nDataPoints; i++)
        {
            if ( y.pointerValue[i] > halfShearHeight ) {
                initialShearFlow.pointerValue[i] = U0*(halfHeight - y.pointerValue[i]);
            } else if ( y.pointerValue[i] < -halfShearHeight ) {
                initialShearFlow.pointerValue[i] = U0*(halfHeight + y.pointerValue[i]);
            } else {
                initialShearFlow.pointerValue[i] = U0*(halfHeight - halfShearHeight/2 - (y.pointerValue[i])*(y.pointerValue[i])/(2*halfShearHeight));
            }
            
            //initialShearFlow.pointerValue[i] = initialShearFlow.pointerValue[i] + 0.001*U0* ( ((GLFloat)rand())/ ((GLFloat) RAND_MAX));
            // Now perturb the shear flow, to get the instability started.
            initialShearFlow.pointerValue[i] = initialShearFlow.pointerValue[i] * (1 + 0.001 * sin( 4*M_PI*x.pointerValue[i]/halfWidth));
        }
        initialShearFlow.differentiationBasis = differentiationBasis;
        GLFunction *zeta = [initialShearFlow differentiateWithOperator: laplacian];
        
        /************************************************************************************************/
        /*		Plop down a float at each grid point													*/
        /************************************************************************************************/
        
        GLFunction *xPosition = [GLFunction functionFromFunction: x];
        GLFunction *yPosition = [GLFunction functionFromFunction: y];
        
        /************************************************************************************************/
        /*		Compute an appropriate time step														*/
        /************************************************************************************************/
        
        CGFloat cfl = 0.25;
        GLFloat timeStep = cfl * xDim.sampleInterval / U0;
        
        GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[zeta, xPosition, yPosition] stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
            GLFunction *psi = [inverseLaplacian transform: yNew[0]];
            GLFunction *advect = [[[psi y] times: [psi differentiateWithOperator: diffJacobianX]] minus: [[psi x] times: [psi differentiateWithOperator: diffJacobianY]]];
            GLFunction *advectSpectral = [advect transformToBasis: @[@(kGLExponentialBasis),@(kGLSineBasis)]];
            GLFunction *f = [advectSpectral plus: [psi differentiateWithOperator: diffLin]];
            
            NSArray *uv = @[[[[psi y] spatialDomain] negate], [[psi x] spatialDomain] ];
            NSArray *xy = @[yNew[1], yNew[2]];
            GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: uv secondOperand: xy];
            
            return @[f, interp.result[0], interp.result[1]];
        }];
        // Tolerances of the float data need to be reduced because they're computed from a linear interpolation of the velocities.
        integrator.absoluteTolerance = @[ @(1e-6), @(1e-3), @(1e-3)];
        
        // Now we create a mutable variable in order to record the evolution of the Gaussian.
        GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [[NSURL fileURLWithPath: [NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) firstObject]] URLByAppendingPathComponent:@"KelvinHelmholtzInstability.nc"] forEquation: equation overwriteExisting: YES];
        
        integrator.shouldDisplayProgress = YES;
        [integrator integrateAlongDimension: tDim withTimeScale: 1.0 file: netcdfFile output: ^(GLScalar *t, NSArray *y) {
            //NSLog(@"Logging day: %f, step size: %f.", (qg.T_QG*rkint.currentTime/86400), rkint.lastStepSize*qg.T_QG);
            
            NSMutableDictionary *scaledVariables = [NSMutableDictionary dictionary];
            
            GLFunction *zeta = [y[0] spatialDomain];
            GLFunction *xPosition = y[1];
            GLFunction *yPosition = y[2];
            GLFunction *psi = [[y[0] differentiateWithOperator: inverseLaplacian] spatialDomain];
            GLFunction *u = [[[psi y] negate] spatialDomain];
            GLFunction *v = [[psi x] spatialDomain];
            
            scaledVariables[@"zeta"] = [zeta scaleVariableBy: 1.0 withUnits: @"1/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"psi"] = [psi scaleVariableBy: 1.0 withUnits: @"m^2/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"u"] = [u scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"v"] = [v scaleVariableBy: 1.0 withUnits: @"m/s" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"x-position"] = [xPosition scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"m"];
            scaledVariables[@"y-position"] = [yPosition scaleVariableBy: 1.0 withUnits: @"m" dimensionsBy: 1.0 units: @"m"];
            
            return scaledVariables;
        }];
        
        NSLog(@"Close the NetCDF file and wrap up");
        
        [equation waitUntilAllOperationsAreFinished];
        
        // The NetCDF file may still be writing data. We need to make sure it finishes before we exit.
        [netcdfFile waitUntilAllOperationsAreFinished];
        [netcdfFile close];
    }
    return 0;
}

