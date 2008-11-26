/*
 * NestedMICA Motif Inference Toolkit
 *
 * Copyright (c) 2004-2007: Genome Research Ltd.
 *
 * NestedMICA is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * or see the on-line version at http://www.gnu.org/copyleft/lgpl.txt
 *
 */

package net.derkholm.nmica.build;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import com.sun.mirror.apt.AnnotationProcessor;
import com.sun.mirror.apt.AnnotationProcessorEnvironment;
import com.sun.mirror.declaration.TypeDeclaration;

public class NMExtraAppProcessor implements AnnotationProcessor {
	private final AnnotationProcessorEnvironment env;
	private File output;
	
	public NMExtraAppProcessor(AnnotationProcessorEnvironment env) {
		this.env = env;
		for (String opt : env.getOptions().keySet()) {
			if (opt.startsWith("-AlaunchpadDir=")) {
				this.output = new File(opt.substring("-AlaunchpadDir=".length()));
			}
		}
	}
	
	public void process() {
		for (TypeDeclaration td : env.getTypeDeclarations()) {
			NMExtraApp nmapp = td.getAnnotation(NMExtraApp.class);
			if (nmapp != null) {
				if (this.output == null) {
					env.getMessager().printError("Missing required option -AlaunchpadDir");
				} else {
					File launchpad = new File(output, nmapp.launchName());
					// env.getMessager().printNotice(String.format("Should generate launchpad %s", launchpad));
					try {
						PrintWriter pw = new PrintWriter(new FileWriter(launchpad));
						pw.printf("#!/bin/sh%n");
						pw.printf("# This file is autogenerated.  If you find problems, please edit NMAppProcessor.java instead%n");
						pw.printf("%n");
						pw.printf("WRAPPER=$0%n");
						pw.printf("APP_DIR=`dirname $WRAPPER`/..%n");
						pw.printf("JVM=%s%n", nmapp.vm() == VirtualMachine.SERVER ? "-server" : "-client");
						pw.printf("MAINCLASS=%s%n", td.getQualifiedName() + "Application");
						pw.printf("%n");
						pw.printf("if ! [ -e \"${APP_DIR}/lib/nmica.jar\" ]; then%n");
						pw.printf("  echo \"Couldn't find nmica.jar.  Have you compiled NestedMICA?\"%n");
						pw.printf("  exit 1%n");
						pw.printf("fi%n");
						pw.printf("%n");
						pw.printf("# If JAVA_HOME is set use that as a JVM pointer.%n");
						pw.printf("%n");
						pw.printf("if [ -z $JAVA_HOME ]; then%n");
						pw.printf("  JAVA_CMD=java%n");
						pw.printf("else%n");
						pw.printf("  JAVA_CMD=${JAVA_HOME}/bin/java%n");
						pw.printf("fi%n");
						pw.printf("%n");
						pw.printf("APP_CLASSPATH=${APP_DIR}/lib/changeless.jar:${APP_DIR}/lib/biojava.jar:${APP_DIR}/lib/bytecode.jar:${APP_DIR}/lib/bjv2-core-0.1.jar:${APP_DIR}/lib/stax-api-1.0.1.jar:${APP_DIR}/lib/wstx-lgpl-3.0.2.jar:${NMICA_DEV_HOME}/lib/nmica.jar:${APP_DIR}/lib/nmica-extra.jar%n");
						pw.printf("%n");
						pw.printf("${JAVA_CMD} ${JVM} ${JVMOPTS} -classpath ${APP_CLASSPATH} -Djava.library.path=${NMICA_HOME}/native -Dchangeless.no_dire_warning=true ${MAINCLASS} \"$@\"%n");

						pw.close();
					} catch (IOException ex) {
						env.getMessager().printError(ex.getMessage());
					}
				}
			}
		}
	}
}
