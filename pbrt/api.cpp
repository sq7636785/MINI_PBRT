
/*
pbrt source code is Copyright(c) 1998-2016
Matt Pharr, Greg Humphreys, and Wenzel Jakob.

This file is part of pbrt.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

// core/api.cpp*
#include "api.h"
#include "parallel.h"
#include "paramset.h"
#include "spectrum.h"
#include "scene.h"
#include "film.h"
#include "medium.h"
#include "stats.h"
#include "sampling.h"

// API Additional Headers
#include "bvh.h"
#include "me/voxel.h"
#include "me/path.h"
#include "me/perspective.h"
#include "me/directlighting.h"
#include "me/point.h"
#include "me/distant.h"
#include "me/curve.h"
#include "me/hair.h"
#include "me/diffuse.h"
#include "me/matte.h"
#include "me/stratified.h"
#include "me/disk.h"
#include "me/sphere.h"
#include "me/infinite.h"
#include "triangle.h"
#include "textures/constant.h"
#include "homogeneous.h"
#include "box.h"
#include "me/random.h"
#include <map>
#include <stdio.h>

namespace pbrt {

	// API Global Variables
	Options PbrtOptions;

	// API Local Classes
	PBRT_CONSTEXPR int MaxTransforms = 2;
	PBRT_CONSTEXPR int StartTransformBits = 1 << 0;
	PBRT_CONSTEXPR int EndTransformBits = 1 << 1;
	PBRT_CONSTEXPR int AllTransformsBits = (1 << MaxTransforms) - 1;
	struct TransformSet {
		// TransformSet Public Methods
		Transform &operator[](int i) {
			CHECK_GE(i, 0);
			CHECK_LT(i, MaxTransforms);
			return t[i];
		}
		const Transform &operator[](int i) const {
			CHECK_GE(i, 0);
			CHECK_LT(i, MaxTransforms);
			return t[i];
		}
		friend TransformSet Inverse(const TransformSet &ts) {
			TransformSet tInv;
			for (int i = 0; i < MaxTransforms; ++i) tInv.t[i] = Inverse(ts.t[i]);
			return tInv;
		}
		bool IsAnimated() const {
			for (int i = 0; i < MaxTransforms - 1; ++i)
				if (t[i] != t[i + 1]) return true;
			return false;
		}

	private:
		Transform t[MaxTransforms];
	};

	struct RenderOptions {
		// RenderOptions Public Methods
		Integrator *MakeIntegrator() const;
		Scene *MakeScene();
		Camera *MakeCamera() const;

		// RenderOptions Public Data
		Float transformStartTime = 0, transformEndTime = 1;
		std::string FilterName = "box";
		ParamSet FilterParams;
		std::string FilmName = "image";
		ParamSet FilmParams;
		std::string SamplerName = "halton";
		ParamSet SamplerParams;
		std::string AcceleratorName = "bvh";
		ParamSet AcceleratorParams;
		std::string IntegratorName = "path";
		ParamSet IntegratorParams;
		std::string CameraName = "perspective";
		ParamSet CameraParams;
		TransformSet CameraToWorld;
		std::map<std::string, std::shared_ptr<Medium>> namedMedia;
		std::vector<std::shared_ptr<Light>> lights;
		std::vector<std::shared_ptr<Primitive>> primitives;
		std::map<std::string, std::vector<std::shared_ptr<Primitive>>> instances;
		std::vector<std::shared_ptr<Primitive>> *currentInstance = nullptr;
		bool haveScatteringMedia = false;
	};

	// MaterialInstance represents both an instance of a material as well as
	// the information required to create another instance of it (possibly with
	// different parameters from the shape).
	struct MaterialInstance {
		MaterialInstance() = default;
		MaterialInstance(const std::string &name, const std::shared_ptr<Material> &mtl,
			ParamSet params)
			: name(name), material(mtl), params(std::move(params)) {
		}

		std::string name;
		std::shared_ptr<Material> material;
		ParamSet params;
	};

	struct GraphicsState {
		// Graphics State Methods
		GraphicsState()
			: floatTextures(std::make_shared<FloatTextureMap>()),
			spectrumTextures(std::make_shared<SpectrumTextureMap>()),
			namedMaterials(std::make_shared<NamedMaterialMap>()) {
			ParamSet empty;
			TextureParams tp(empty, empty, *floatTextures, *spectrumTextures);
			std::shared_ptr<Material> mtl(CreateMatteMaterial(tp));
			currentMaterial = std::make_shared<MaterialInstance>("matte", mtl, ParamSet());
		}
		std::shared_ptr<Material> GetMaterialForShape(const ParamSet &geomParams);
		MediumInterface CreateMediumInterface();

		// Graphics State
		std::string currentInsideMedium, currentOutsideMedium;

		// Updated after book publication: floatTextures, spectrumTextures, and
		// namedMaterials are all implemented using a "copy on write" approach
		// for more efficient GraphicsState management.  When state is pushed
		// in pbrtAttributeBegin(), we don't immediately make a copy of these
		// maps, but instead record that each one is shared.  Only if an item
		// is added to one is a unique copy actually made.
		using FloatTextureMap = std::map<std::string, std::shared_ptr<Texture<Float>>>;
		std::shared_ptr<FloatTextureMap> floatTextures;
		bool floatTexturesShared = false;

		using SpectrumTextureMap = std::map<std::string, std::shared_ptr<Texture<Spectrum>>>;
		std::shared_ptr<SpectrumTextureMap> spectrumTextures;
		bool spectrumTexturesShared = false;

		using NamedMaterialMap = std::map<std::string, std::shared_ptr<MaterialInstance>>;
		std::shared_ptr<NamedMaterialMap> namedMaterials;
		bool namedMaterialsShared = false;

		std::shared_ptr<MaterialInstance> currentMaterial;
		ParamSet areaLightParams;
		std::string areaLight;
		bool reverseOrientation = false;
	};

	STAT_MEMORY_COUNTER("Memory/TransformCache", transformCacheBytes);
	STAT_PERCENT("Scene/TransformCache hits", nTransformCacheHits, nTransformCacheLookups);
	STAT_INT_DISTRIBUTION("Scene/Probes per TransformCache lookup", transformCacheProbes);

	// Note: TransformCache has been reimplemented and has a slightly different
	// interface compared to the version described in the third edition of
	// Physically Based Rendering.  The new version is more efficient in both
	// space and memory, which is helpful for highly complex scenes.
	//
	// The new implementation uses a hash table to store Transforms (rather
	// than a std::map, which generally uses a red-black tree).  Further,
	// it doesn't always store the inverse of the transform; if a caller
	// wants the inverse as well, they are responsible for storing it.
	//
	// The hash table size is always a power of two, allowing for the use of a
	// bitwise AND to turn hash values into table offsets.  Quadratic probing
	// is used when there is a hash collision.
	class TransformCache {
	public:
		TransformCache()
			: hashTable(512), hashTableOccupancy(0) {
		}

		// TransformCache Public Methods
		Transform *Lookup(const Transform &t) {
			++nTransformCacheLookups;

			int offset = Hash(t) & (hashTable.size() - 1);
			int step = 1;
			while (true) {
				// Keep looking until we find the Transform or determine that
				// it's not present.
				if (!hashTable[offset] || *hashTable[offset] == t)
					break;
				// Advance using quadratic probing.
				offset = (offset + step * step) & (hashTable.size() - 1);
				++step;
			}
			ReportValue(transformCacheProbes, step);
			Transform *tCached = hashTable[offset];
			if (tCached)
				++nTransformCacheHits;
			else {
				tCached = arena.Alloc<Transform>();
				*tCached = t;
				Insert(tCached);
			}
			return tCached;
		}

		void Clear() {
			transformCacheBytes += arena.TotalAllocated() + hashTable.size() * sizeof(Transform *);
			hashTable.resize(512);
			hashTable.clear();
			hashTableOccupancy = 0;
			arena.Reset();
		}

	private:
		void Insert(Transform *tNew);
		void Grow();

		static uint64_t Hash(const Transform &t) {
			const char *ptr = (const char *)(&t.GetMatrix());
			size_t size = sizeof(Matrix4x4);
			uint64_t hash = 14695981039346656037ull;
			while (size > 0) {
				hash ^= *ptr;
				hash *= 1099511628211ull;
				++ptr;
				--size;
			}
			return hash;
		}

		// TransformCache Private Data
		std::vector<Transform *> hashTable;
		int hashTableOccupancy;
		MemoryArena arena;
	};

	void TransformCache::Insert(Transform *tNew) {
		if (++hashTableOccupancy == hashTable.size() / 2)
			Grow();

		int offset = Hash(*tNew) & (hashTable.size() - 1);
		int step = 1;
		while (true) {
			if (hashTable[offset] == nullptr) {
				hashTable[offset] = tNew;
				return;
			}
			// Advance using quadratic probing.
			offset = (offset + step * step) & (hashTable.size() - 1);
			++step;
		}
	}

	void TransformCache::Grow() {
		std::vector<Transform *> newTable(2 * hashTable.size());
		LOG(INFO) << "Growing transform cache hash table to " << newTable.size();

		// Insert current elements into newTable.
		for (Transform *tEntry : hashTable) {
			if (!tEntry) continue;

			int offset = Hash(*tEntry) & (newTable.size() - 1);
			int step = 1;
			while (true) {
				if (newTable[offset] == nullptr) {
					newTable[offset] = tEntry;
					break;
				}
				// Advance using quadratic probing.
				offset = (offset + step * step) & (hashTable.size() - 1);
				++step;
			}
		}

		std::swap(hashTable, newTable);
	}


	// API Static Data
	enum class APIState { Uninitialized, OptionsBlock, WorldBlock };
	static APIState currentApiState = APIState::Uninitialized;
	static TransformSet curTransform;
	static uint32_t activeTransformBits = AllTransformsBits;
	static std::map<std::string, TransformSet> namedCoordinateSystems;
	static std::unique_ptr<RenderOptions> renderOptions;
	static GraphicsState graphicsState;
	static std::vector<GraphicsState> pushedGraphicsStates;
	static std::vector<TransformSet> pushedTransforms;
	static std::vector<uint32_t> pushedActiveTransformBits;
	static TransformCache transformCache;
	int catIndentCount = 0;

	// API Forward Declarations
	std::vector<std::shared_ptr<Shape>> MakeShapes(const std::string &name,
		const Transform *ObjectToWorld,
		const Transform *WorldToObject,
		bool reverseOrientation,
		const ParamSet &paramSet);

	// API Macros
#define VERIFY_INITIALIZED(func)                           \
    if (!(PbrtOptions.cat || PbrtOptions.toPly) &&           \
        currentApiState == APIState::Uninitialized) {        \
        Error(                                             \
            "pbrtInit() must be before calling \"%s()\". " \
            "Ignoring.",                                   \
            func);                                         \
        return;                                            \
    } else /* swallow trailing semicolon */
#define VERIFY_OPTIONS(func)                             \
    VERIFY_INITIALIZED(func);                            \
    if (!(PbrtOptions.cat || PbrtOptions.toPly) &&       \
        currentApiState == APIState::WorldBlock) {       \
        Error(                                           \
            "Options cannot be set inside world block; " \
            "\"%s\" not allowed.  Ignoring.",            \
            func);                                       \
        return;                                          \
    } else /* swallow trailing semicolon */
#define VERIFY_WORLD(func)                                   \
    VERIFY_INITIALIZED(func);                                \
    if (!(PbrtOptions.cat || PbrtOptions.toPly) &&           \
        currentApiState == APIState::OptionsBlock) {         \
        Error(                                               \
            "Scene description must be inside world block; " \
            "\"%s\" not allowed. Ignoring.",                 \
            func);                                           \
        return;                                              \
    } else /* swallow trailing semicolon */
#define FOR_ACTIVE_TRANSFORMS(expr)           \
    for (int i = 0; i < MaxTransforms; ++i)   \
        if (activeTransformBits & (1 << i)) { \
            expr                              \
        }
#define WARN_IF_ANIMATED_TRANSFORM(func)                             \
    do {                                                             \
        if (curTransform.IsAnimated())                               \
            Warning(                                                 \
                "Animated transformations set; ignoring for \"%s\" " \
                "and using the start transform only",                \
                func);                                               \
    } while (false) /* swallow trailing semicolon */

	// Object Creation Function Definitions
	std::vector<std::shared_ptr<Shape>> MakeShapes(const std::string &name,
		const Transform *object2world,
		const Transform *world2object,
		bool reverseOrientation,
		const ParamSet &paramSet) {
		std::vector<std::shared_ptr<Shape>> shapes;
		std::shared_ptr<Shape> s;
		if (name == "sphere")
			s = CreateSphereShape(object2world, world2object, reverseOrientation,
				paramSet);
		// Create remaining single _Shape_ types
		else if (name == "disk")
			s = CreateDiskShape(object2world, world2object, reverseOrientation,
				paramSet);
		if (s != nullptr) shapes.push_back(s);
		else if (name == "curve") {
			shapes = CreateCurveShape(object2world, world2object,
				reverseOrientation, paramSet);
		}
		else if (name == "trianglemesh") {
			if (PbrtOptions.toPly) {
				int nvi;
				const int *vi = paramSet.FindInt("indices", &nvi);

				if (nvi < 500) {
					// It's a small mesh; don't bother with a PLY file after all.
					printf("%*sShape \"%s\" ", catIndentCount, "", name.c_str());
					paramSet.Print(catIndentCount);
					printf("\n");
				} else {
					static int count = 1;
					const char *plyPrefix =
						getenv("PLY_PREFIX") ? getenv("PLY_PREFIX") : "mesh";
					std::string fn = StringPrintf("%s_%05d.ply", plyPrefix, count++);

					int npi, nuvi, nsi, nni;
					const Point3f *P = paramSet.FindPoint3f("P", &npi);
					const Point2f *uvs = paramSet.FindPoint2f("uv", &nuvi);
					if (!uvs) uvs = paramSet.FindPoint2f("st", &nuvi);
					std::vector<Point2f> tempUVs;
					if (!uvs) {
						const Float *fuv = paramSet.FindFloat("uv", &nuvi);
						if (!fuv) fuv = paramSet.FindFloat("st", &nuvi);
						if (fuv) {
							nuvi /= 2;
							tempUVs.reserve(nuvi);
							for (int i = 0; i < nuvi; ++i)
								tempUVs.push_back(Point2f(fuv[2 * i], fuv[2 * i + 1]));
							uvs = &tempUVs[0];
						}
					}
					const Normal3f *N = paramSet.FindNormal3f("N", &nni);
					const Vector3f *S = paramSet.FindVector3f("S", &nsi);
					int nfi;
					const int *faceIndices = paramSet.FindInt("faceIndices", &nfi);
					if (faceIndices) CHECK_EQ(nfi, nvi / 3);


					ParamSet ps = paramSet;
					ps.EraseInt("indices");
					ps.ErasePoint3f("P");
					ps.ErasePoint2f("uv");
					ps.ErasePoint2f("st");
					ps.EraseNormal3f("N");
					ps.EraseVector3f("S");
					ps.EraseInt("faceIndices");

					printf("%*sShape \"plymesh\" \"string filename\" \"%s\" ",
						catIndentCount, "", fn.c_str());
					ps.Print(catIndentCount);
					printf("\n");
				}
			} else
				shapes = CreateTriangleMeshShape(object2world, world2object,
					reverseOrientation, paramSet,
					&*graphicsState.floatTextures);
		}
		else
			Warning("Shape \"%s\" unknown.", name.c_str());
		return shapes;
	}

	STAT_COUNTER("Scene/Materials created", nMaterialsCreated);

	std::shared_ptr<Material> MakeMaterial(const std::string &name,
		const TextureParams &mp) {
		Material *material = nullptr;
		if (name == "" || name == "none")
			return nullptr;
		else if (name == "matte")
			material = CreateMatteMaterial(mp);
		else if (name == "hair")
			material = CreateHairMaterial(mp);
		else {
			Warning("Material \"%s\" unknown. Using \"matte\".", name.c_str());
			material = CreateMatteMaterial(mp);
		}

		if ((name == "subsurface" || name == "kdsubsurface") &&
			(renderOptions->IntegratorName != "path" &&
			(renderOptions->IntegratorName != "volpath")))
			Warning(
				"Subsurface scattering material \"%s\" used, but \"%s\" "
				"integrator doesn't support subsurface scattering. "
				"Use \"path\" or \"volpath\".",
				name.c_str(), renderOptions->IntegratorName.c_str());

		mp.ReportUnused();
		if (!material) Error("Unable to create material \"%s\"", name.c_str());
		else ++nMaterialsCreated;
		return std::shared_ptr<Material>(material);
	}

	std::shared_ptr<Texture<Float>> MakeFloatTexture(const std::string &name,
		const Transform &tex2world,
		const TextureParams &tp) {
		Texture<Float> *tex = nullptr;
		if (name == "constant")
			tex = CreateConstantFloatTexture(tex2world, tp);
		else
			Warning("Float texture \"%s\" unknown.", name.c_str());
		tp.ReportUnused();
		return std::shared_ptr<Texture<Float>>(tex);
	}

	std::shared_ptr<Texture<Spectrum>> MakeSpectrumTexture(
		const std::string &name, const Transform &tex2world,
		const TextureParams &tp) {
		Texture<Spectrum> *tex = nullptr;
		if (name == "constant")
			tex = CreateConstantSpectrumTexture(tex2world, tp);
		else
			Warning("Spectrum texture \"%s\" unknown.", name.c_str());
		tp.ReportUnused();
		return std::shared_ptr<Texture<Spectrum>>(tex);
	}

	std::shared_ptr<Medium> MakeMedium(const std::string &name,
		const ParamSet &paramSet,
		const Transform &medium2world) {
		Float sig_a_rgb[3] = { .0011f, .0024f, .014f },
			sig_s_rgb[3] = { 2.55f, 3.21f, 3.77f };
		Spectrum sig_a = Spectrum::FromRGB(sig_a_rgb),
			sig_s = Spectrum::FromRGB(sig_s_rgb);
		std::string preset = paramSet.FindOneString("preset", "");
		bool found = GetMediumScatteringProperties(preset, &sig_a, &sig_s);
		if (preset != "" && !found)
			Warning("Material preset \"%s\" not found.  Using defaults.",
				preset.c_str());
		Float scale = paramSet.FindOneFloat("scale", 1.f);
		Float g = paramSet.FindOneFloat("g", 0.0f);
		sig_a = paramSet.FindOneSpectrum("sigma_a", sig_a) * scale;
		sig_s = paramSet.FindOneSpectrum("sigma_s", sig_s) * scale;
		Medium *m = NULL;
		if (name == "homogeneous") {
			m = new HomogeneousMedium(sig_a, sig_s, g);
		} else
			Warning("Medium \"%s\" unknown.", name.c_str());
		paramSet.ReportUnused();
		return std::shared_ptr<Medium>(m);
	}

	std::shared_ptr<Light> MakeLight(const std::string &name,
		const ParamSet &paramSet,
		const Transform &light2world,
		const MediumInterface &mediumInterface) {
		std::shared_ptr<Light> light;
		if (name == "point")
			light =
			CreatePointLight(light2world, mediumInterface.outside, paramSet);
		else if (name == "distant")
			light = CreateDistantLight(light2world, paramSet);
		else if (name == "infinite" || name == "exinfinite")
			light = CreateInfiniteLight(light2world, paramSet);
		else
			Warning("Light \"%s\" unknown.", name.c_str());
		paramSet.ReportUnused();
		return light;
	}

	std::shared_ptr<AreaLight> MakeAreaLight(const std::string &name,
		const Transform &light2world,
		const MediumInterface &mediumInterface,
		const ParamSet &paramSet,
		const std::shared_ptr<Shape> &shape) {
		std::shared_ptr<AreaLight> area;
		if (name == "area" || name == "diffuse")
			area = CreateDiffuseAreaLight(light2world, mediumInterface.outside,
				paramSet, shape);
		else
			Warning("Area light \"%s\" unknown.", name.c_str());
		paramSet.ReportUnused();
		return area;
	}

	std::shared_ptr<Primitive> MakeAccelerator(
		const std::string &name,
		std::vector<std::shared_ptr<Primitive>> prims,
		const ParamSet &paramSet) {
		std::shared_ptr<Primitive> accel;
		if (name == "bvh")
			accel = CreateBVHAccelerator(std::move(prims), paramSet);
		else
			Warning("Accelerator \"%s\" unknown.", name.c_str());
		paramSet.ReportUnused();
		return accel;
	}

	Camera *MakeCamera(const std::string &name, const ParamSet &paramSet,
		const TransformSet &cam2worldSet, Float transformStart,
		Float transformEnd, Film *film) {
		Camera *camera = nullptr;
		MediumInterface mediumInterface = graphicsState.CreateMediumInterface();
		static_assert(MaxTransforms == 2,
			"TransformCache assumes only two transforms");
		Transform *cam2world[2] = {
			transformCache.Lookup(cam2worldSet[0]),
			transformCache.Lookup(cam2worldSet[1])
		};
		AnimatedTransform animatedCam2World(cam2world[0], transformStart,
			cam2world[1], transformEnd);
		if (name == "perspective")
			camera = CreatePerspectiveCamera(paramSet, animatedCam2World, film,
				mediumInterface.outside);
		else
			Warning("Camera \"%s\" unknown.", name.c_str());
		paramSet.ReportUnused();
		return camera;
	}

	std::shared_ptr<Sampler> MakeSampler(const std::string &name,
		const ParamSet &paramSet,
		const Film *film) {
		Sampler *sampler = nullptr;
		if (name == "stratified")
			sampler = CreateStratifiedSampler(paramSet);
		else if (name == "random")
			sampler = CreateRandomSampler(paramSet);
		else
			Warning("Sampler \"%s\" unknown.", name.c_str());
		paramSet.ReportUnused();
		return std::shared_ptr<Sampler>(sampler);
	}

	std::unique_ptr<Filter> MakeFilter(const std::string &name,
		const ParamSet &paramSet) {
		Filter *filter = nullptr;

		if (name == "box")
			filter = CreateBoxFilter(paramSet);
		else {
			Error("Filter \"%s\" unknown.", name.c_str());
			exit(1);
		}
		paramSet.ReportUnused();
		return std::unique_ptr<Filter>(filter);
	}

	Film *MakeFilm(const std::string &name, const ParamSet &paramSet,
		std::unique_ptr<Filter> filter) {
		Film *film = nullptr;
		if (name == "image")
			film = CreateFilm(paramSet, std::move(filter));
		else
			Warning("Film \"%s\" unknown.", name.c_str());
		paramSet.ReportUnused();
		return film;
	}

	// API Function Definitions
	void pbrtInit(const Options &opt) {
		PbrtOptions = opt;
		// API Initialization
		if (currentApiState != APIState::Uninitialized)
			Error("pbrtInit() has already been called.");
		currentApiState = APIState::OptionsBlock;
		renderOptions.reset(new RenderOptions);
		graphicsState = GraphicsState();
		catIndentCount = 0;

		// General \pbrt Initialization
		SampledSpectrum::Init();
		ParallelInit();  // Threads must be launched before the profiler is
						 // initialized.
		InitProfiler();
	}

	void pbrtCleanup() {
		// API Cleanup
		if (currentApiState == APIState::Uninitialized)
			Error("pbrtCleanup() called without pbrtInit().");
		else if (currentApiState == APIState::WorldBlock)
			Error("pbrtCleanup() called while inside world block.");
		currentApiState = APIState::Uninitialized;
		ParallelCleanup();
		CleanupProfiler();
	}

	void pbrtIdentity() {
		VERIFY_INITIALIZED("Identity");
		FOR_ACTIVE_TRANSFORMS(curTransform[i] = Transform();)
			if (PbrtOptions.cat || PbrtOptions.toPly)
				printf("%*sIdentity\n", catIndentCount, "");
	}

	void pbrtTranslate(Float dx, Float dy, Float dz) {
		VERIFY_INITIALIZED("Translate");
		FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] *
			Translate(Vector3f(dx, dy, dz));)
			if (PbrtOptions.cat || PbrtOptions.toPly)
				printf("%*sTranslate %.9g %.9g %.9g\n", catIndentCount, "", dx, dy,
					dz);
	}

	void pbrtTransform(Float tr[16]) {
		VERIFY_INITIALIZED("Transform");
		FOR_ACTIVE_TRANSFORMS(
			curTransform[i] = Transform(Matrix4x4(
				tr[0], tr[4], tr[8], tr[12], tr[1], tr[5], tr[9], tr[13], tr[2],
				tr[6], tr[10], tr[14], tr[3], tr[7], tr[11], tr[15]));)
			if (PbrtOptions.cat || PbrtOptions.toPly) {
				printf("%*sTransform [ ", catIndentCount, "");
				for (int i = 0; i < 16; ++i) printf("%.9g ", tr[i]);
				printf(" ]\n");
			}
	}

	void pbrtConcatTransform(Float tr[16]) {
		VERIFY_INITIALIZED("ConcatTransform");
		FOR_ACTIVE_TRANSFORMS(
			curTransform[i] =
			curTransform[i] *
			Transform(Matrix4x4(tr[0], tr[4], tr[8], tr[12], tr[1], tr[5],
				tr[9], tr[13], tr[2], tr[6], tr[10], tr[14],
				tr[3], tr[7], tr[11], tr[15]));)
			if (PbrtOptions.cat || PbrtOptions.toPly) {
				printf("%*sConcatTransform [ ", catIndentCount, "");
				for (int i = 0; i < 16; ++i) printf("%.9g ", tr[i]);
				printf(" ]\n");
			}
	}

	void pbrtRotate(Float angle, Float dx, Float dy, Float dz) {
		VERIFY_INITIALIZED("Rotate");
		FOR_ACTIVE_TRANSFORMS(curTransform[i] =
			curTransform[i] *
			Rotate(angle, Vector3f(dx, dy, dz));)
			if (PbrtOptions.cat || PbrtOptions.toPly)
				printf("%*sRotate %.9g %.9g %.9g %.9g\n", catIndentCount, "", angle,
					dx, dy, dz);
	}

	void pbrtScale(Float sx, Float sy, Float sz) {
		VERIFY_INITIALIZED("Scale");
		FOR_ACTIVE_TRANSFORMS(curTransform[i] =
			curTransform[i] * Scale(sx, sy, sz);)
			if (PbrtOptions.cat || PbrtOptions.toPly)
				printf("%*sScale %.9g %.9g %.9g\n", catIndentCount, "", sx, sy, sz);
	}

	void pbrtLookAt(Float ex, Float ey, Float ez, Float lx, Float ly, Float lz,
		Float ux, Float uy, Float uz) {
		VERIFY_INITIALIZED("LookAt");
		Transform lookAt =
			LookAt(Point3f(ex, ey, ez), Point3f(lx, ly, lz), Vector3f(ux, uy, uz));
		FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * lookAt;);
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf(
				"%*sLookAt %.9g %.9g %.9g\n%*s%.9g %.9g %.9g\n"
				"%*s%.9g %.9g %.9g\n",
				catIndentCount, "", ex, ey, ez, catIndentCount + 8, "", lx, ly, lz,
				catIndentCount + 8, "", ux, uy, uz);
	}

	void pbrtCoordinateSystem(const std::string &name) {
		VERIFY_INITIALIZED("CoordinateSystem");
		namedCoordinateSystems[name] = curTransform;
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sCoordinateSystem \"%s\"\n", catIndentCount, "",
				name.c_str());
	}

	void pbrtCoordSysTransform(const std::string &name) {
		VERIFY_INITIALIZED("CoordSysTransform");
		if (namedCoordinateSystems.find(name) != namedCoordinateSystems.end())
			curTransform = namedCoordinateSystems[name];
		else
			Warning("Couldn't find named coordinate system \"%s\"", name.c_str());
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sCoordSysTransform \"%s\"\n", catIndentCount, "",
				name.c_str());
	}

	void pbrtActiveTransformAll() {
		activeTransformBits = AllTransformsBits;
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sActiveTransform All\n", catIndentCount, "");
	}

	void pbrtActiveTransformEndTime() {
		activeTransformBits = EndTransformBits;
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sActiveTransform EndTime\n", catIndentCount, "");
	}

	void pbrtActiveTransformStartTime() {
		activeTransformBits = StartTransformBits;
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sActiveTransform StartTime\n", catIndentCount, "");
	}

	void pbrtTransformTimes(Float start, Float end) {
		VERIFY_OPTIONS("TransformTimes");
		renderOptions->transformStartTime = start;
		renderOptions->transformEndTime = end;
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sTransformTimes %.9g %.9g\n", catIndentCount, "", start,
				end);
	}

	void pbrtPixelFilter(const std::string &name, const ParamSet &params) {
		VERIFY_OPTIONS("PixelFilter");
		renderOptions->FilterName = name;
		renderOptions->FilterParams = params;
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sPixelFilter \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtFilm(const std::string &type, const ParamSet &params) {
		VERIFY_OPTIONS("Film");
		renderOptions->FilmParams = params;
		renderOptions->FilmName = type;
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sFilm \"%s\" ", catIndentCount, "", type.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtSampler(const std::string &name, const ParamSet &params) {
		VERIFY_OPTIONS("Sampler");
		renderOptions->SamplerName = name;
		renderOptions->SamplerParams = params;
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sSampler \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtAccelerator(const std::string &name, const ParamSet &params) {
		VERIFY_OPTIONS("Accelerator");
		renderOptions->AcceleratorName = name;
		renderOptions->AcceleratorParams = params;
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sAccelerator \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtIntegrator(const std::string &name, const ParamSet &params) {
		VERIFY_OPTIONS("Integrator");
		renderOptions->IntegratorName = name;
		renderOptions->IntegratorParams = params;
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sIntegrator \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtCamera(const std::string &name, const ParamSet &params) {
		VERIFY_OPTIONS("Camera");
		renderOptions->CameraName = name;
		renderOptions->CameraParams = params;
		renderOptions->CameraToWorld = Inverse(curTransform);
		namedCoordinateSystems["camera"] = renderOptions->CameraToWorld;
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sCamera \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtMakeNamedMedium(const std::string &name, const ParamSet &params) {
		VERIFY_INITIALIZED("MakeNamedMedium");
		WARN_IF_ANIMATED_TRANSFORM("MakeNamedMedium");
		std::string type = params.FindOneString("type", "");
		if (type == "")
			Error("No parameter string \"type\" found in MakeNamedMedium");
		else {
			std::shared_ptr<Medium> medium =
				MakeMedium(type, params, curTransform[0]);
			if (medium) renderOptions->namedMedia[name] = medium;
		}
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sMakeNamedMedium \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtMediumInterface(const std::string &insideName,
		const std::string &outsideName) {
		VERIFY_INITIALIZED("MediumInterface");
		graphicsState.currentInsideMedium = insideName;
		graphicsState.currentOutsideMedium = outsideName;
		renderOptions->haveScatteringMedia = true;
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sMediumInterface \"%s\" \"%s\"\n", catIndentCount, "",
				insideName.c_str(), outsideName.c_str());
	}

	void pbrtWorldBegin() {
		VERIFY_OPTIONS("WorldBegin");
		currentApiState = APIState::WorldBlock;
		for (int i = 0; i < MaxTransforms; ++i) curTransform[i] = Transform();
		activeTransformBits = AllTransformsBits;
		namedCoordinateSystems["world"] = curTransform;
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("\n\nWorldBegin\n\n");
	}

	void pbrtAttributeBegin() {
		VERIFY_WORLD("AttributeBegin");
		pushedGraphicsStates.push_back(graphicsState);
		graphicsState.floatTexturesShared = graphicsState.spectrumTexturesShared =
			graphicsState.namedMaterialsShared = true;
		pushedTransforms.push_back(curTransform);
		pushedActiveTransformBits.push_back(activeTransformBits);
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("\n%*sAttributeBegin\n", catIndentCount, "");
			catIndentCount += 4;
		}
	}

	void pbrtAttributeEnd() {
		VERIFY_WORLD("AttributeEnd");
		if (!pushedGraphicsStates.size()) {
			Error(
				"Unmatched pbrtAttributeEnd() encountered. "
				"Ignoring it.");
			return;
		}
		graphicsState = std::move(pushedGraphicsStates.back());
		pushedGraphicsStates.pop_back();
		curTransform = pushedTransforms.back();
		pushedTransforms.pop_back();
		activeTransformBits = pushedActiveTransformBits.back();
		pushedActiveTransformBits.pop_back();
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			catIndentCount -= 4;
			printf("%*sAttributeEnd\n", catIndentCount, "");
		}
	}

	void pbrtTransformBegin() {
		VERIFY_WORLD("TransformBegin");
		pushedTransforms.push_back(curTransform);
		pushedActiveTransformBits.push_back(activeTransformBits);
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sTransformBegin\n", catIndentCount, "");
			catIndentCount += 4;
		}
	}

	void pbrtTransformEnd() {
		VERIFY_WORLD("TransformEnd");
		if (!pushedTransforms.size()) {
			Error(
				"Unmatched pbrtTransformEnd() encountered. "
				"Ignoring it.");
			return;
		}
		curTransform = pushedTransforms.back();
		pushedTransforms.pop_back();
		activeTransformBits = pushedActiveTransformBits.back();
		pushedActiveTransformBits.pop_back();
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			catIndentCount -= 4;
			printf("%*sTransformEnd\n", catIndentCount, "");
		}
	}

	void pbrtTexture(const std::string &name, const std::string &type,
		const std::string &texname, const ParamSet &params) {
		VERIFY_WORLD("Texture");
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sTexture \"%s\" \"%s\" \"%s\" ", catIndentCount, "",
				name.c_str(), type.c_str(), texname.c_str());
			params.Print(catIndentCount);
			printf("\n");
			return;
		}

		TextureParams tp(params, params, *graphicsState.floatTextures,
			*graphicsState.spectrumTextures);
		if (type == "float") {
			// Create _Float_ texture and store in _floatTextures_
			if (graphicsState.floatTextures->find(name) !=
				graphicsState.floatTextures->end())
				Warning("Texture \"%s\" being redefined", name.c_str());
			WARN_IF_ANIMATED_TRANSFORM("Texture");
			std::shared_ptr<Texture<Float>> ft =
				MakeFloatTexture(texname, curTransform[0], tp);
			if (ft) {
				// TODO: move this to be a GraphicsState method, also don't
				// provide direct floatTextures access?
				if (graphicsState.floatTexturesShared) {
					graphicsState.floatTextures =
						std::make_shared<GraphicsState::FloatTextureMap>(*graphicsState.floatTextures);
					graphicsState.floatTexturesShared = false;
				}
				(*graphicsState.floatTextures)[name] = ft;
			}
		} else if (type == "color" || type == "spectrum") {
			// Create _color_ texture and store in _spectrumTextures_
			if (graphicsState.spectrumTextures->find(name) !=
				graphicsState.spectrumTextures->end())
				Warning("Texture \"%s\" being redefined", name.c_str());
			WARN_IF_ANIMATED_TRANSFORM("Texture");
			std::shared_ptr<Texture<Spectrum>> st =
				MakeSpectrumTexture(texname, curTransform[0], tp);
			if (st) {
				if (graphicsState.spectrumTexturesShared) {
					graphicsState.spectrumTextures =
						std::make_shared<GraphicsState::SpectrumTextureMap>(*graphicsState.spectrumTextures);
					graphicsState.spectrumTexturesShared = false;
				}
				(*graphicsState.spectrumTextures)[name] = st;
			}
		} else
			Error("Texture type \"%s\" unknown.", type.c_str());
	}

	void pbrtMaterial(const std::string &name, const ParamSet &params) {
		VERIFY_WORLD("Material");
		ParamSet emptyParams;
		TextureParams mp(params, emptyParams, *graphicsState.floatTextures,
			*graphicsState.spectrumTextures);
		std::shared_ptr<Material> mtl = MakeMaterial(name, mp);
		graphicsState.currentMaterial =
			std::make_shared<MaterialInstance>(name, mtl, params);

		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sMaterial \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtMakeNamedMaterial(const std::string &name, const ParamSet &params) {
		VERIFY_WORLD("MakeNamedMaterial");
		// error checking, warning if replace, what to use for transform?
		ParamSet emptyParams;
		TextureParams mp(params, emptyParams, *graphicsState.floatTextures,
			*graphicsState.spectrumTextures);
		std::string matName = mp.FindString("type");
		WARN_IF_ANIMATED_TRANSFORM("MakeNamedMaterial");
		if (matName == "")
			Error("No parameter string \"type\" found in MakeNamedMaterial");

		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sMakeNamedMaterial \"%s\" ", catIndentCount, "",
				name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		} else {
			std::shared_ptr<Material> mtl = MakeMaterial(matName, mp);
			if (graphicsState.namedMaterials->find(name) !=
				graphicsState.namedMaterials->end())
				Warning("Named material \"%s\" redefined.", name.c_str());
			if (graphicsState.namedMaterialsShared) {
				graphicsState.namedMaterials =
					std::make_shared<GraphicsState::NamedMaterialMap>(*graphicsState.namedMaterials);
				graphicsState.namedMaterialsShared = false;
			}
			(*graphicsState.namedMaterials)[name] =
				std::make_shared<MaterialInstance>(matName, mtl, params);
		}
	}

	void pbrtNamedMaterial(const std::string &name) {
		VERIFY_WORLD("NamedMaterial");
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sNamedMaterial \"%s\"\n", catIndentCount, "", name.c_str());
			return;
		}

		auto iter = graphicsState.namedMaterials->find(name);
		if (iter == graphicsState.namedMaterials->end()) {
			Error("NamedMaterial \"%s\" unknown.", name.c_str());
			return;
		}
		graphicsState.currentMaterial = iter->second;
	}

	void pbrtLightSource(const std::string &name, const ParamSet &params) {
		VERIFY_WORLD("LightSource");
		WARN_IF_ANIMATED_TRANSFORM("LightSource");
		MediumInterface mi = graphicsState.CreateMediumInterface();
		std::shared_ptr<Light> lt = MakeLight(name, params, curTransform[0], mi);
		if (!lt)
			Error("LightSource: light type \"%s\" unknown.", name.c_str());
		else
			renderOptions->lights.push_back(lt);
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sLightSource \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtAreaLightSource(const std::string &name, const ParamSet &params) {
		VERIFY_WORLD("AreaLightSource");
		graphicsState.areaLight = name;
		graphicsState.areaLightParams = params;
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sAreaLightSource \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}
	}

	void pbrtShape(const std::string &name, const ParamSet &params) {
		VERIFY_WORLD("Shape");
		std::vector<std::shared_ptr<Primitive>> prims;
		std::vector<std::shared_ptr<AreaLight>> areaLights;
		if (PbrtOptions.cat || (PbrtOptions.toPly && name != "trianglemesh")) {
			printf("%*sShape \"%s\" ", catIndentCount, "", name.c_str());
			params.Print(catIndentCount);
			printf("\n");
		}

		if (!curTransform.IsAnimated()) {
			// Initialize _prims_ and _areaLights_ for static shape

			// Create shapes for shape _name_
			Transform *ObjToWorld = transformCache.Lookup(curTransform[0]);
			Transform *WorldToObj = transformCache.Lookup(Inverse(curTransform[0]));
			std::vector<std::shared_ptr<Shape>> shapes =
				MakeShapes(name, ObjToWorld, WorldToObj,
					graphicsState.reverseOrientation, params);
			if (shapes.empty()) return;
			std::shared_ptr<Material> mtl = graphicsState.GetMaterialForShape(params);
			params.ReportUnused();
			MediumInterface mi = graphicsState.CreateMediumInterface();
			prims.reserve(shapes.size());
			for (auto s : shapes) {
				// Possibly create area light for shape
				std::shared_ptr<AreaLight> area;
				if (graphicsState.areaLight != "") {
					area = MakeAreaLight(graphicsState.areaLight, curTransform[0],
						mi, graphicsState.areaLightParams, s);
					if (area) areaLights.push_back(area);
				}
				prims.push_back(
					std::make_shared<GeometricPrimitive>(s, mtl, area, mi));
			}
		} else {
			// Initialize _prims_ and _areaLights_ for animated shape

			// Create initial shape or shapes for animated shape
			if (graphicsState.areaLight != "")
				Warning(
					"Ignoring currently set area light when creating "
					"animated shape");
			Transform *identity = transformCache.Lookup(Transform());
			std::vector<std::shared_ptr<Shape>> shapes = MakeShapes(
				name, identity, identity, graphicsState.reverseOrientation, params);
			if (shapes.empty()) return;

			// Create _GeometricPrimitive_(s) for animated shape
			std::shared_ptr<Material> mtl = graphicsState.GetMaterialForShape(params);
			params.ReportUnused();
			MediumInterface mi = graphicsState.CreateMediumInterface();
			prims.reserve(shapes.size());
			for (auto s : shapes)
				prims.push_back(
					std::make_shared<GeometricPrimitive>(s, mtl, nullptr, mi));

			// Create single _TransformedPrimitive_ for _prims_

			// Get _animatedObjectToWorld_ transform for shape
			static_assert(MaxTransforms == 2,
				"TransformCache assumes only two transforms");
			Transform *ObjToWorld[2] = {
				transformCache.Lookup(curTransform[0]),
				transformCache.Lookup(curTransform[1])
			};
			AnimatedTransform animatedObjectToWorld(
				ObjToWorld[0], renderOptions->transformStartTime, ObjToWorld[1],
				renderOptions->transformEndTime);
			if (prims.size() > 1) {
				std::shared_ptr<Primitive> bvh = std::make_shared<BVHAccel>(prims);
				prims.clear();
				prims.push_back(bvh);
			}
			prims[0] = std::make_shared<TransformedPrimitive>(
				prims[0], animatedObjectToWorld);
		}
		// Add _prims_ and _areaLights_ to scene or current instance
		if (renderOptions->currentInstance) {
			if (areaLights.size())
				Warning("Area lights not supported with object instancing");
			renderOptions->currentInstance->insert(
				renderOptions->currentInstance->end(), prims.begin(), prims.end());
		} else {
			renderOptions->primitives.insert(renderOptions->primitives.end(),
				prims.begin(), prims.end());
			if (areaLights.size())
				renderOptions->lights.insert(renderOptions->lights.end(),
					areaLights.begin(), areaLights.end());
		}
	}

	// Attempt to determine if the ParamSet for a shape may provide a value for
	// its material's parameters. Unfortunately, materials don't provide an
	// explicit representation of their parameters that we can query and
	// cross-reference with the parameter values available from the shape.
	//
	// Therefore, we'll apply some "heuristics".
	bool shapeMaySetMaterialParameters(const ParamSet &ps) {
		for (const auto &param : ps.textures)
			// Any texture other than one for an alpha mask is almost certainly
			// for a Material (or is unused!).
			if (param->name != "alpha" && param->name != "shadowalpha")
				return true;

		// Special case spheres, which are the most common non-mesh primitive.
		for (const auto &param : ps.floats)
			if (param->nValues == 1 && param->name != "radius")
				return true;

		// Extra special case strings, since plymesh uses "filename", curve "type",
		// and loopsubdiv "scheme".
		for (const auto &param : ps.strings)
			if (param->nValues == 1 && param->name != "filename" &&
				param->name != "type" && param->name != "scheme")
				return true;

		// For all other parameter types, if there is a single value of the
		// parameter, assume it may be for the material. This should be valid
		// (if conservative), since no materials currently take array
		// parameters.
		for (const auto &param : ps.bools)
			if (param->nValues == 1)
				return true;
		for (const auto &param : ps.ints)
			if (param->nValues == 1)
				return true;
		for (const auto &param : ps.point2fs)
			if (param->nValues == 1)
				return true;
		for (const auto &param : ps.vector2fs)
			if (param->nValues == 1)
				return true;
		for (const auto &param : ps.point3fs)
			if (param->nValues == 1)
				return true;
		for (const auto &param : ps.vector3fs)
			if (param->nValues == 1)
				return true;
		for (const auto &param : ps.normals)
			if (param->nValues == 1)
				return true;
		for (const auto &param : ps.spectra)
			if (param->nValues == 1)
				return true;

		return false;
	}

	std::shared_ptr<Material> GraphicsState::GetMaterialForShape(
		const ParamSet &shapeParams) {
		CHECK(currentMaterial);
		if (shapeMaySetMaterialParameters(shapeParams)) {
			// Only create a unique material for the shape if the shape's
			// parameters are (apparently) going to provide values for some of
			// the material parameters.
			TextureParams mp(shapeParams, currentMaterial->params, *floatTextures,
				*spectrumTextures);
			return MakeMaterial(currentMaterial->name, mp);
		} else
			return currentMaterial->material;
	}

	MediumInterface GraphicsState::CreateMediumInterface() {
		MediumInterface m;
		if (currentInsideMedium != "") {
			if (renderOptions->namedMedia.find(currentInsideMedium) !=
				renderOptions->namedMedia.end())
				m.inside = renderOptions->namedMedia[currentInsideMedium].get();
			else
				Error("Named medium \"%s\" undefined.",
					currentInsideMedium.c_str());
		}
		if (currentOutsideMedium != "") {
			if (renderOptions->namedMedia.find(currentOutsideMedium) !=
				renderOptions->namedMedia.end())
				m.outside = renderOptions->namedMedia[currentOutsideMedium].get();
			else
				Error("Named medium \"%s\" undefined.",
					currentOutsideMedium.c_str());
		}
		return m;
	}

	void pbrtReverseOrientation() {
		VERIFY_WORLD("ReverseOrientation");
		graphicsState.reverseOrientation = !graphicsState.reverseOrientation;
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sReverseOrientation\n", catIndentCount, "");
	}

	void pbrtObjectBegin(const std::string &name) {
		VERIFY_WORLD("ObjectBegin");
		pbrtAttributeBegin();
		if (renderOptions->currentInstance)
			Error("ObjectBegin called inside of instance definition");
		renderOptions->instances[name] = std::vector<std::shared_ptr<Primitive>>();
		renderOptions->currentInstance = &renderOptions->instances[name];
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sObjectBegin \"%s\"\n", catIndentCount, "", name.c_str());
	}

	STAT_COUNTER("Scene/Object instances created", nObjectInstancesCreated);

	void pbrtObjectEnd() {
		VERIFY_WORLD("ObjectEnd");
		if (!renderOptions->currentInstance)
			Error("ObjectEnd called outside of instance definition");
		renderOptions->currentInstance = nullptr;
		pbrtAttributeEnd();
		++nObjectInstancesCreated;
		if (PbrtOptions.cat || PbrtOptions.toPly)
			printf("%*sObjectEnd\n", catIndentCount, "");
	}

	STAT_COUNTER("Scene/Object instances used", nObjectInstancesUsed);

	void pbrtObjectInstance(const std::string &name) {
		VERIFY_WORLD("ObjectInstance");
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sObjectInstance \"%s\"\n", catIndentCount, "", name.c_str());
			return;
		}

		// Perform object instance error checking
		if (renderOptions->currentInstance) {
			Error("ObjectInstance can't be called inside instance definition");
			return;
		}
		if (renderOptions->instances.find(name) == renderOptions->instances.end()) {
			Error("Unable to find instance named \"%s\"", name.c_str());
			return;
		}
		std::vector<std::shared_ptr<Primitive>> &in =
			renderOptions->instances[name];
		if (in.empty()) return;
		++nObjectInstancesUsed;
		if (in.size() > 1) {
			// Create aggregate for instance _Primitive_s
			std::shared_ptr<Primitive> accel(
				MakeAccelerator(renderOptions->AcceleratorName, std::move(in),
					renderOptions->AcceleratorParams));
			if (!accel) accel = std::make_shared<BVHAccel>(in);
			in.clear();
			in.push_back(accel);
		}
		static_assert(MaxTransforms == 2,
			"TransformCache assumes only two transforms");
		// Create _animatedInstanceToWorld_ transform for instance
		Transform *InstanceToWorld[2] = {
			transformCache.Lookup(curTransform[0]),
			transformCache.Lookup(curTransform[1])
		};
		AnimatedTransform animatedInstanceToWorld(
			InstanceToWorld[0], renderOptions->transformStartTime,
			InstanceToWorld[1], renderOptions->transformEndTime);
		std::shared_ptr<Primitive> prim(
			std::make_shared<TransformedPrimitive>(in[0], animatedInstanceToWorld));
		renderOptions->primitives.push_back(prim);
	}

	void pbrtWorldEnd() {
		VERIFY_WORLD("WorldEnd");
		// Ensure there are no pushed graphics states
		while (pushedGraphicsStates.size()) {
			Warning("Missing end to pbrtAttributeBegin()");
			pushedGraphicsStates.pop_back();
			pushedTransforms.pop_back();
		}
		while (pushedTransforms.size()) {
			Warning("Missing end to pbrtTransformBegin()");
			pushedTransforms.pop_back();
		}

		// Create scene and render
		if (PbrtOptions.cat || PbrtOptions.toPly) {
			printf("%*sWorldEnd\n", catIndentCount, "");
		} else {
			std::unique_ptr<Integrator> integrator(renderOptions->MakeIntegrator());
			std::unique_ptr<Scene> scene(renderOptions->MakeScene());

			// This is kind of ugly; we directly override the current profiler
			// state to switch from parsing/scene construction related stuff to
			// rendering stuff and then switch it back below. The underlying
			// issue is that all the rest of the profiling system assumes
			// hierarchical inheritance of profiling state; this is the only
			// place where that isn't the case.
			CHECK_EQ(CurrentProfilerState(), ProfToBits(Prof::SceneConstruction));
			ProfilerState = ProfToBits(Prof::IntegratorRender);

			if (scene && integrator) integrator->Render(*scene);

			CHECK_EQ(CurrentProfilerState(), ProfToBits(Prof::IntegratorRender));
			ProfilerState = ProfToBits(Prof::SceneConstruction);
		}

		// Clean up after rendering. Do this before reporting stats so that
		// destructors can run and update stats as needed.
		graphicsState = GraphicsState();
		transformCache.Clear();
		currentApiState = APIState::OptionsBlock;
		
		renderOptions.reset(new RenderOptions);

		if (!PbrtOptions.cat && !PbrtOptions.toPly) {
			MergeWorkerThreadStats();
			ReportThreadStats();
			if (!PbrtOptions.quiet) {
				PrintStats(stdout);
				ReportProfilerResults(stdout);
				ClearStats();
				ClearProfiler();
			}
		}

		for (int i = 0; i < MaxTransforms; ++i) curTransform[i] = Transform();
		activeTransformBits = AllTransformsBits;
		namedCoordinateSystems.erase(namedCoordinateSystems.begin(),
			namedCoordinateSystems.end());
	}

#define  VOXEL_INFO

	Scene *RenderOptions::MakeScene() {
		auto shapeCopy = primitives;
		std::shared_ptr<Primitive> accelerator =
			MakeAccelerator(AcceleratorName, std::move(primitives), AcceleratorParams);
		if (!accelerator) accelerator = std::make_shared<BVHAccel>(primitives);

		std::shared_ptr<Volume> volume = std::make_shared<Volume>(accelerator->WorldBound(), 50.0, 15, 64 * 64);
		clock_t s = clock();
		volume->ConstructVolume();
  		volume->CalculateVoxel(shapeCopy, true);
		clock_t e = clock();
		std::cout << "voxelization time: " << (e - s) << " ms" << std::endl;

#ifdef VOXEL_INFO

		std::cout << "volume voxel minDim " << volume->GetMinDimDelta() << " " << volume->GetMaxDimDelta() << std::endl;

#endif

		Scene *scene = new Scene(accelerator, volume, lights);

		bool isMultiple = false;

		if (isMultiple) {
			s = clock();
			scene->InitHairBSDF(true);

			scene->VolumeIndirectLight(3000000);

			//scene->VolumeIrrandiance();
			//scene->VolumeSHRadiance();

			e = clock();
			//std::cout << "irrandiance time: " << (e - s) << " ms" << std::endl;
			std::cout << "volume randiance time: " << (e - s) << " ms" << std::endl;

			s = clock();
			scene->volume->BoxFilterSHC();
			e = clock();
			//std::cout << "irrandiance time: " << (e - s) << " ms" << std::endl;
			std::cout << "box filter time: " << (e - s) << " ms" << std::endl;

			//bsdfMatrix
			s = clock();
			//scene->VolumeBSDFMatrix(40000, 256);
			scene->volume->LoadBsdfMatrix("./bsdf_matrix/15_40000_256.txt");
			e = clock();
			std::cout << "volume bsdf Matrix time " << (e - s) << " ms" << std::endl;
		}
		//		volume->LoadData("volumeFull_10-10-10.txt");
		//volume->SaveData("radiance_50-50-50Irrandiance.txt");

		/************************************************************************/
		/*  test volume                                                         */
		/************************************************************************/
// 		bool t1 = volume->ConstructVolume();
// 		bool t2 = volume->CalculateVoxel(shapeCopy);
// 		Point3f startPos(1.3, 2.1, -1.5);
// 		Point3f endPos(4.2, 3.5, 0.5);
// 		Spectrum tr = volume->Tr(startPos, endPos);
		/************************************************************************/
		/*                                                                      */
		/************************************************************************/


		// Erase primitives and lights from _RenderOptions_
		primitives.clear();
		lights.clear();
		return scene;
	}

	Integrator *RenderOptions::MakeIntegrator() const {
		std::shared_ptr<const Camera> camera(MakeCamera());
		if (!camera) {
			Error("Unable to create camera");
			return nullptr;
		}

		std::shared_ptr<Sampler> sampler =
			MakeSampler(SamplerName, SamplerParams, camera->film);
		if (!sampler) {
			Error("Unable to create sampler.");
			return nullptr;
		}

		Integrator *integrator = nullptr;
		if (IntegratorName == "directlighting")
			integrator =
			CreateDirectLightingIntegrator(IntegratorParams, sampler, camera);
		else if (IntegratorName == "path")
			integrator =
			CreatePathIntegrator(IntegratorParams, sampler, camera);
		else {
			Error("Integrator \"%s\" unknown.", IntegratorName.c_str());
			return nullptr;
		}

		if (renderOptions->haveScatteringMedia && IntegratorName != "volpath" &&
			IntegratorName != "bdpt" && IntegratorName != "mlt") {
			Warning(
				"Scene has scattering media but \"%s\" integrator doesn't support "
				"volume scattering. Consider using \"volpath\", \"bdpt\", or "
				"\"mlt\".", IntegratorName.c_str());
		}

		IntegratorParams.ReportUnused();
		// Warn if no light sources are defined
		if (lights.empty())
			Warning(
				"No light sources defined in scene; "
				"rendering a black image.");
		return integrator;
	}

	Camera *RenderOptions::MakeCamera() const {
		std::unique_ptr<Filter> filter = MakeFilter(FilterName, FilterParams);
		Film *film = MakeFilm(FilmName, FilmParams, std::move(filter));
		if (!film) {
			Error("Unable to create film.");
			return nullptr;
		}
		Camera *camera = pbrt::MakeCamera(CameraName, CameraParams, CameraToWorld,
			renderOptions->transformStartTime,
			renderOptions->transformEndTime, film);
		return camera;
	}

	/**/

	void TestHairF() {
		RNG rng;
		Vector3f wo = UniformSampleSphere({ rng.UniformFloat(), rng.UniformFloat() });
		for (Float beta_m = .1; beta_m < 1; beta_m += .2) {
			for (Float beta_n = .1; beta_n < 1; beta_n += .2) {
				// Estimate reflected uniform incident radiance from hair
				Spectrum sum = 0.f;
				int count = 300000;
				for (int i = 0; i < count; ++i) {
					Float h = -1 + 2. * rng.UniformFloat();
					Spectrum sigma_a = 0.f;
					HairBSDF hair(h, 1.55, sigma_a, beta_m, beta_n, 0.f);
					Vector3f wi = UniformSampleSphere(
						{ rng.UniformFloat(), rng.UniformFloat() });
					sum += hair.f(wo, wi) * AbsCosTheta(wi);
				}
				Float avg = sum.y() / (count * UniformSpherePdf());
				std::cout << "f: " << avg << std::endl;
			}
		}
	}

	void TestHairSampleF() {
		RNG rng;
		Vector3f wo = UniformSampleSphere({ rng.UniformFloat(), rng.UniformFloat() });
		for (Float beta_m = .1; beta_m < 1; beta_m += .2) {
			for (Float beta_n = .1; beta_n < 1; beta_n += .2) {
				Spectrum sum = 0.f;
				int count = 300000;
				for (int i = 0; i < count; ++i) {
					Float h = -1 + 2. * rng.UniformFloat();
					Spectrum sigma_a = 0.f;
					HairBSDF hair(h, 1.55, sigma_a, beta_m, beta_n, 0.f);

					Vector3f wi;
					Float pdf;
					Point2f u = { rng.UniformFloat(), rng.UniformFloat() };
					Spectrum f = hair.Sample_f(wo, &wi, u, &pdf, nullptr);
					if (pdf > 0) sum += f * AbsCosTheta(wi) / pdf;
				}
				Float avg = sum.y() / count;
				std::cout << "sample f: " << avg << std::endl;
			}
		}
	}

	void TestHairPdf() {
		RNG rng;
		for (Float beta_m = .1; beta_m < 1; beta_m += .2)
			for (Float beta_n = .4; beta_n < 1; beta_n += .2) {
				int count = 10000;
				for (int i = 0; i < count; ++i) {
					// Check _HairBSDF::Sample\_f()_ sample weight
					Float h = -1 + 2 * rng.UniformFloat();
					Spectrum sigma_a = 0;
					HairBSDF hair(h, 1.55, sigma_a, beta_m, beta_n, 0.f);
					Vector3f wo = UniformSampleSphere(
						{ rng.UniformFloat(), rng.UniformFloat() });
					Vector3f wi;
					Float pdf;
					Point2f u = { rng.UniformFloat(), rng.UniformFloat() };
					Spectrum f = hair.Sample_f(wo, &wi, u, &pdf, nullptr);
					if (pdf > 0) {
						// Verify that hair BSDF sample weight is close to 1 for
						// _wi_
						f = (f * AbsCosTheta(wi) / pdf);
						std::cout << "pdf weight: " << f.y() << std::endl;
					}
				}
			}
	}

	

}  // namespace pbrt
