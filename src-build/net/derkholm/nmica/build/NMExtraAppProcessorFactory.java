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

import java.util.Arrays;
import java.util.Collection;
import java.util.Set;

import com.sun.mirror.apt.AnnotationProcessor;
import com.sun.mirror.apt.AnnotationProcessorEnvironment;
import com.sun.mirror.apt.AnnotationProcessorFactory;
import com.sun.mirror.declaration.AnnotationTypeDeclaration;

public class MXplorAppProcessorFactory implements AnnotationProcessorFactory {
	public AnnotationProcessor getProcessorFor(
			Set<AnnotationTypeDeclaration> decl, 
			AnnotationProcessorEnvironment env) 
	{
		return new MXplorAppProcessor(env);
	}

	public Collection<String> supportedAnnotationTypes() {
		return Arrays.asList(
			new String[] {
				MXplorApp.class.getName(), 
				MXplorAppProcessor.class.getName()});
	}

	public Collection<String> supportedOptions() {
		return Arrays.asList("-AlaunchpadDir");
	}
}
