package net.derkholm.nmica.extra.app.seq.nextgen;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.bjv2.util.cli.App;

@NMExtraApp(launchName = "ngcount", vm = VirtualMachine.SERVER)
@App(overview = "Output sequencing depth inside a window.", generateStub = true)
public class CountReads {

}
